from selenium import webdriver
from selenium.webdriver.common.by import By
import time
from openai import OpenAI
import pandas as pd
import tqdm
import os
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail

browser = webdriver.Chrome()
client = OpenAI()

#Scrapes NSF SBIR table until specified year, default is 2019
def scrape_table(data, year):
    print("Scraping table...")
    tbody = browser.find_element(By.TAG_NAME, "tbody")
    for row in tqdm.tqdm(tbody.find_elements(By.TAG_NAME, "tr")):
        row_class = row.get_attribute("class")
        if row_class == "even" or row_class == "odd":
            columns = row.find_elements(By.TAG_NAME, "td")

            if columns[5].text.split("/")[2] == year:
                return data
            else:
                data["Link"].append(columns[0].find_element(By.TAG_NAME, "a").get_attribute("href"))
                data["Company"].append(columns[0].text)
                data["Location"].append(columns[1].text)
                data["Title"].append(columns[3].text)
                data["Date"].append(columns[5].text)

                #the last entry in the table
                if(columns[2] == "1058288"):
                    return data
    return data

#Runs the table scraper
def run_scrape_table(time_sleep=1.5, year_to_stop="2022", filename="nsf.csv"):
    data = {
        "Link": [],
        "Company":[],
        "Location":[],
        "Title": [],
        "Date":[],
    }
    browser.get("https://seedfund.nsf.gov/awardees/history/?yr=All&p1&p2&kw&pp=-1&gv&st")
    time.sleep(time_sleep)

    data = scrape_table(data, year_to_stop)

    df = pd.DataFrame(data)
    df.to_csv(filename)

def filter_by_state(states=[", MD", ", VA"], filename="nsf.csv"):
    usecols = ['Link','Company','Location','Title','Date']
    df = pd.read_csv(filename, usecols=usecols)
    mask = df['Location'].str.contains('|'.join(states))
    local_companies = df[mask]
    local_companies.to_csv(f"data/local_relevant.csv")

def filter_by_relevance(system_message, prompt_example, assistant_example, filename="nsf.csv"):
    # GOES IN MAIN:
    #filter_by_relevance(
    #     system_message="Determine if the following project is very directly relevant to Computer Science. Simply output '1' if it does and '0' if it doesn't.",
    #     prompt_example="Repurposing RNA processing factors to edit RNA for the treatment of Dravet syndrome",
    #     assistant_example="0",
    #     filename=filename
    # )
    usecols = ['Link','Company','Location','Title','Date']
    df = pd.read_csv(filename, usecols=usecols)

    relevant_rows = []

    print("Filtering by relevance...")
    for index, row in tqdm.tqdm(df.iterrows()):
        if (index == 500):
            break
        prompt = row['Title']
        print(prompt)
        response = client.chat.completions.create(
            model="gpt-4-0125-preview",
            messages=[
                {"role": "system", "content": system_message},
                {"role": "user", "content": prompt_example},
                {"role": "assistant", "content": assistant_example},
                {"role": "user", "content": prompt},
            ],
        )
        ai_response = response.choices[0].message.content
        print(ai_response)
        if '1' in ai_response:
            relevant_rows.append(row)
    
    new_df = pd.DataFrame(relevant_rows)
    new_df.to_csv("data/airelevant.csv")


def prepare_email(row, system_message):
    browser.get(row["Link"])
    time.sleep(3)
    try:
        abstract = browser.find_elements(By.CLASS_NAME, 'abstractText')[1].text.replace('\n', ' ')
        address = browser.find_elements(By.TAG_NAME, 'a')[18].get_attribute('href')
        summarize_response = client.chat.completions.create(
                    model="gpt-3.5-turbo-0125",
                    messages=[
                        {"role": "system", "content": "Please summarize this paragraph into 3 sentences while maintaining as much original intent and meaning as possible. In your response, identify detailed information about project goals, technical detail, and, if possible, the concrete outcome the project aims to achieve."},
                        {"role": "user", "content": abstract},
                    ],
                )
        summarized = summarize_response.choices[0].message.content
        email_response = client.chat.completions.create(
                model="gpt-4-0125-preview",
                messages=[
                    {"role": "system", "content": system_message},
                    {"role": "user", "content": summarized},
                ],
            )
        email = email_response.choices[0].message.content
        return address, email
    except Exception as e:
        print(f"Erorr processing {row['Link']}:  {e}\n")
        return 0, 0
    
def format_name(name):
    response = client.chat.completions.create(
                    model="gpt-3.5-turbo-0125",
                    messages=[
                        {"role": "system", "content": "Please format the following into title case and remove any extensions like LLC, Inc, etc."},
                        {"role": "user", "content": name},
                    ],
                )
    return response.choices[0].message.content

def send_emails():
    df = pd.read_csv("data/rest.csv")
    for index, row in tqdm.tqdm(df.iterrows()):
        name = format_name(row["Company"])
        subject = f"Request: Seeking Internship Opportunity at {name}"
        sender_email = "hudson@hudsontao.com"

        email, expression_of_interest = prepare_email(
            row,
            system_message="You are a job application expert helping me write a cold email. I will give you a summary of the company I'm applying to, and your task is generate a single super brief sentence expressing my interest in their work. The goal is to show that I have a general idea of what their company does. Do not surround your output with quotations. Keep the language in layman's terms and be very general. The sentence should have no more than 1 comma separator.",
        )
        if (email != 0):

            remote_version = f"""
                    <div>
                        <p>
                            Hello,
                        </p>
                    </div>
                    <div>
                        <br>
                    </div>
                    <div>
                        <p>
                            I'm Hudson Tao, an undergrad CS major at the University of Maryland. I'm writing to ask for an internship opportunity at {name} this summer. {expression_of_interest} I have previous work experience as a researcher and course creator for educational AI startup Learn Prompting, as well as a software developer for startup LifeSct LLC, a high-end biotech product supplier. I would really like to dedicate my work and effort to something that really matters, and I believe that an internship opportunity at your company would let me achieve that. <a href='hudsontao.com'>Here is a link to my resume</a>, which goes into more details about my technical background. Thank you so much for your consideration.
                        </p>
                    </div>
                    <div>
                        <br>
                    </div>
                    <div>
                        <p>
                            Very respectfully,
                        </p>
                    </div>
                    <div>
                        <br>
                    </div>
                    <div>
                        <p>
                            Hudson Tao
                        </p>
                    </div>
                    <div>
                        <p>
                            240-743-0735
                        </p>
                    </div>
                    <div>
                        <br>
                    </div>
                    <div>
                        <p>
                            Please note: Given the considerable distance, I am looking to work remotely.
                        </p>
                    </div>
                """
            local_version = f"""
                    <div>
                        <p>
                            Hello,
                        </p>
                    </div>
                    <div>
                        <br>
                    </div>
                    <div>
                        <p>
                            I'm Hudson Tao, an undergrad CS major at the University of Maryland. I'm writing to ask for an internship opportunity at {name} this summer. {expression_of_interest} I have previous work experience as a researcher and course creator for educational AI startup Learn Prompting, as well as a software developer for startup LifeSct LLC, a high-end biotech product supplier. I would really like to dedicate my work and effort to something that really matters, and I believe that an internship opportunity at your company would let me achieve that. <a href='hudsontao.com'>Here is a link to my resume</a>, which goes into more details about my technical background. Thank you so much for your consideration.
                        </p>
                    </div>
                    <div>
                        <br>
                    </div>
                    <div>
                        <p>
                            Very respectfully,
                        </p>
                    </div>
                    <div>
                        <p>
                            Hudson Tao
                        </p>
                    </div>
                    <div>
                        <p>
                            240-743-0735
                        </p>
                    </div>
                    <div>
                        <br>
                    </div>
                    <div>
                        <p>
                            Please note: I am willing to work both remotely and in-person.
                        </p>
                    </div>
                """
            local = ["POTOMAC", "SILVER SPRING", "ROCKVILLE", "LANHAM", "GAITHERSBURG", "TAKOMA PARK"]
            location = row["Location"]
            sendee_email = email
            if any(word in location for word in local):
                message = Mail(
                    from_email=sender_email,
                    to_emails=sendee_email,
                    subject=subject,
                    html_content=local_version
                )
            else:
                message = Mail(
                    from_email=sender_email,
                    to_emails=sendee_email,
                    subject=subject,
                    html_content=remote_version
                )
            try:
                sg = SendGridAPIClient(os.environ.get('SENDGRID_API_KEY'))
                sg.send(message)
                print(f"Email sent to {email}\n")
            except Exception as e:
                print(f"Error sending: {e}\n")
        else:
            print(f"Email {index} not sent.\n")


if __name__ == "__main__":
    filename = "data/nsf.csv"

    send_emails()

    browser.close()