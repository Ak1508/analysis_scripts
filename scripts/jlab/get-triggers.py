#!/usr/local/bin/python3.7

import json
import requests
import re
import datetime
import numpy as np
import matplotlib.pyplot as plt

elogs     = requests.get('https://logbooks.jlab.org/api/elog/entries?book=HCLOG&title=End_of_Run&field=lognumber&field=title&field=body&field=created&limit=75')
jelogs    = elogs.json()
elEntries = jelogs['data']['entries']

elTitle = []; elTimeStamp = []; elBody = []; elContent = []; dateTime = []
runNumber = []; numEvents = []; liveTime = []

for index, entry in enumerate(elEntries) :
    elTitle.append(elEntries[index]['title'])
    elTimeStamp.append(elEntries[index]['created'])
    elBody.append(elEntries[index]['body'])
    elContent.append(elBody[index]['content'])
    dateTime.append(datetime.datetime(int (elTimeStamp[index]['string'][:4]), int (elTimeStamp[index]['string'][5:7]), int (elTimeStamp[index]['string'][8:10]), int(elTimeStamp[index]['string'][11:13]), int(elTimeStamp[index]['string'][14:16])))
    runNumber.append(int (re.findall(r'\d+', elTitle[index])[0]))
    numEvents.append(int (re.findall(r'\d+', elContent[index].partition('Event number = ')[2])[0]))
    liveTime.append(float (re.findall(r'\d+\.\d+', elContent[index].partition('Live Time: ')[2])[0]))

plt.figure('Run Number vs. Number of Events')
plt.plot(runNumber, numEvents, 'md')
plt.xlabel('Run Number')
plt.ylabel('Number of Events')
plt.show()

plt.figure('Run Number vs. Live Time')
plt.plot(runNumber, liveTime, 'bo')
plt.xlabel('Run Number')
plt.ylabel('Live Time (%)')
plt.show()
