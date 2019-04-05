#!/usr/bin/python

import json
import requests
import re
import datetime

elogs  = requests.get('https://logbooks.jlab.org/api/elog/entries?book=HCLOG&title=End_of_Run&field=lognumber&field=title&field=body&field=created&limit=1')
jelogs = elogs.json()

elEntries   = jelogs['data']['entries'][0]
elTitle     = elEntries['title']
elTimeStamp = elEntries['created']
elBody      = elEntries['body']
elContent   = elBody['content']

dateTime = datetime.datetime(int (elTimeStamp['string'][:4]), int (elTimeStamp['string'][5:7]), int (elTimeStamp['string'][8:10]), int(elTimeStamp['string'][11:13]), int(elTimeStamp['string'][14:16]))

runNumber = int (re.findall(r'\d+', elTitle)[0])
numEvents = int (re.findall(r'\d+', elContent.partition('Event number = ')[2])[0])





