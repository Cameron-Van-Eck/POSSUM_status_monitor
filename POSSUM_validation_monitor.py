#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 09:03:40 2024

Automatic POSSUM validation status checker.
Scrapes the validation status from Emil's script on the ATNF servers,
then parses the return for how old the SBs are.
Sends the status to the #validation Slack channel, with an alert if there are
any particularly old SBs.

@author: cvaneck
"""

import argparse
import requests
import re
import os
import pyvo as vo
from astropy.time import Time
import datetime

def parse_script():
    """This function grabs the output of Emil's validation status script, and
    parses it to get the ages of all SBs awaiting validation.
    Returns a dictionary with keys = ages and values = SB numbers. 
    
    Deprecated since Emil's tool is no longer accessible."""
    
    r = requests.get('https://www.atnf.csiro.au/research/RACS/ASKAPview/status.py?proj=AS203')
    lines = r.text.split('\n')
    sb_lines = [ x for x in lines if 'DEPOSITED' in x]
    deposit_ages = [float(re.findall('DEPOSITED (\d+\.\d*) days ago',x)[0]) for x in sb_lines ]

    sb_numbers = [re.findall('(\d+) not yet released',x)[0] for x in sb_lines ]
    
    current_sbs = dict(zip(sb_numbers,deposit_ages))

    
    
    return current_sbs


def get_unvalidated_SBs():
    session=vo.tap.TAPService('https://casda.csiro.au/casda_vo_tools/tap')
    result=session.search("SELECT sbid,event_date,event_type,obs_program,project_code FROM casda.observation_event WHERE project_code = 'AS203' ORDER BY event_date").to_table()
    
    sbid_status = {}
    sbid_time = {}
    for item in result:
        sbid = int(item['sbid'])
        event = item['event_type']
        sbid_status[sbid] = event
        sbid_time[sbid] = item["event_date"]
    
    sbids = list(sbid_status.keys())
    sbids.sort()
    
    unreleased = []
    ages = {}
    t_now = Time(datetime.datetime.now(), scale='utc')
    
    for sbid in sbids:
        if sbid_status[sbid] not in ["UPDATED", "RELEASED","REJECTED", "DELETED", "VALIDATED"]:
            unreleased.append(sbid)
            event_time = Time(sbid_time[sbid].replace("T", " ").replace("Z", ""))
            ages[sbid] = (t_now - event_time).sec / 60.0 / 60.0 / 24.0
    
    
    return ages




def notify_slack(webhook,current_sbs,state_dict):

    #This is the age at which a warning is sent to the channel about potential lockouts.
    warning_age_threshold= 24     

    if len(current_sbs) == 0:
        message=('This is a biweekly automatic update on validation status.\n'
                 'No SBs have been found awaiting validation.\n'
                 "Either we're working very efficiently, no new SBs have been deposited"
                 " or something in Cameron's script has broken. Someone should double check on CASDA.")
        
    else:

        EMU_sb = [x for x in state_dict if (state_dict[x][0] == 'EMU')  ]
        EMU_ages = [current_sbs[x] for x in EMU_sb]

        WALLABY_sb = [x for x in state_dict if (state_dict[x][0] == 'WALLABY')  ]
        WALLABY_ages = [current_sbs[x] for x in WALLABY_sb]


        sb_count=len(current_sbs)

        message=('This is a biweekly automatic update on validation status.\n'
                f'There are currently {sb_count} SBs awaiting validation.\n')


        if len(EMU_sb) == 0:
            message+='There are no EMU SBs awaiting validation.\n'
            oldest_age_EMU=0
        else:
            oldest_age_EMU = max(EMU_ages)
            message+= f'The oldest EMU SB is {oldest_age_EMU:.1f} days old.\n'
            
        if len(WALLABY_sb) == 0:
            message+='There are no WALLABY SBs awaiting validation.\n'
        else:
            oldest_age_WALLABY = max(WALLABY_ages)
            message+= f'The oldest WALLABY SB is {oldest_age_WALLABY:.1f} days old.\n'
                


        #Send warning if SBs are close to deadline
        if oldest_age_EMU > warning_age_threshold:
            message+=(f'\n @channel *Warning: at least one EMU SB is older than {warning_age_threshold} days.* '
                      'Please validate before it hits 30 days to prevent an observation lockout.\n'
                      f'EMU SBs above warning threshold: {[x for x in EMU_sb if current_sbs[x] > warning_age_threshold ]}\n')
    
        if len(EMU_sb) > 0:
            message+='\nEMU SBs:\n'
        
                #Sort by states.
            sb_waiting = [x for x in state_dict if (state_dict[x][1] == 'WAITING') and (state_dict[x][0] == 'EMU')  ]
            sb_ready = [x for x in state_dict if (state_dict[x][1] == 'VALIDATED')  and (state_dict[x][0] == 'EMU')   ]
            sb_reject = [x for x in state_dict if (state_dict[x][1] == 'REJECTED')  and (state_dict[x][0] == 'EMU')  ]
        
            if len(sb_ready) > 0:
                message+=f'    Ready for validation: {sb_ready}\n'
            if len(sb_reject) > 0:
                message+=f'    Rejected by EMU: {sb_reject}\n'
            if len(sb_waiting) > 0:
                message+=f'    Waiting for EMU validation: {sb_waiting}\n'
            

        if len(WALLABY_sb) > 0:
            message+='\nWALLABY SBs:\n'
        
                #Sort by states.
            sb_WALLABY = [x for x in state_dict if (state_dict[x][0] == 'WALLABY')  ]
    #        sb_ready = [x for x in state_dict if (state_dict[x][1] == 'VALIDATED')  and (state_dict[x][0] == 'WALLABY')   ]
    #        sb_reject = [x for x in state_dict if (state_dict[x][1] == 'REJECTED')  and (state_dict[x][0] == 'WALLABY')  ]
        
            if len(sb_WALLABY) > 0:
                message+=f'    Ready for validation: {sb_WALLABY}\n'
            # if len(sb_reject) > 0:
            #     message+=f'    Rejected by WALLABY: {sb_reject}\n'
            # if len(sb_waiting) > 0:
            #     message+=f'    Waiting for WALLABY validation: {sb_waiting}\n'
        
    
    requests.post(webhook,json={"parse":"full",'text':message})
    
    
    
    
def check_validation_state(sb_list):
    
    session=vo.tap.TAPService('https://casda.csiro.au/casda_vo_tools/tap')

    state_dict={}
    for sb in sb_list:
        query=("SELECT TOP 100 "
               "sbid,event_date,event_type,obs_program,project_code "
               "FROM casda.observation_event "
               f"WHERE sbid = '{sb}' AND (project_code = 'AS202' OR project_code = 'AS201')")
        result=session.search(query).to_table()
        if len(result) == 0:
            raise Exception(f"SB {sb} doesn't exist in CASDA?")

        if result[0]['project_code'] == 'AS201':
            survey='EMU'
        elif result[0]['project_code'] == 'AS202':
            survey='WALLABY'
        else:
            raise Exception('SB {sb} not associated with EMU or WALLABY???')
        
        if 'REJECTED' in result['event_type']:
            state_dict[sb] = (survey, ' REJECTED')
        elif 'RELEASED' in result['event_type']:
            state_dict[sb] = (survey, 'VALIDATED')
        else:
            state_dict[sb] = (survey, 'WAITING')

    return state_dict
    

def run():
    
    descStr = """
    Check the POSSUM validation status. Requires a text file with a Slack 
    webhook in order to notify Slack.
    """

    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("webhook_file", metavar="webhook.txt",
                        help="Text file with Slack webhook URL")
    args = parser.parse_args()

    if not os.path.exists(args.webhook_file):
        raise Exception(f'Cannot find webhook file at {args.webhook_file}')
    
    with open(args.webhook_file,'r') as f:
        webhook=f.read()
    webhook=webhook.strip()
    
    
    current_sbs=get_unvalidated_SBs()
    state_dict=check_validation_state(list(current_sbs))
    notify_slack(webhook,current_sbs,state_dict)







if __name__ == "__main__":
    run()
