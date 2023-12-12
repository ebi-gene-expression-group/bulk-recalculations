import sys
import logging
import os
import re

# log handler script for Snakefile-sorting-hat

l = logging.getLogger(__name__)
l.setLevel(logging.DEBUG)
fh = logging.FileHandler(os.environ.get('LOG_PATH'), mode='w')
l.addHandler(fh)

jobs = {}

def log_handler(msg):

    if "jobid" in msg:
        if msg['jobid'] in jobs:
            if 'level' in msg:
                jobs[msg['jobid']]['status'] = msg['level']
        else:
            job = { 'name': msg['name'], 'status': 'created' }
            if 'accession' in msg['wildcards']:
                job['accession'] = msg['wildcards']['accession']
            if 'log' in msg and len(msg['log']) > 0:
                job['log'] = msg['log'][0]
            jobs[msg['jobid']] = job

        for k, v in msg.items():
            print(f"Key: {k}")
            print(f"Value: {v}")

        jobs_to_remove = []
        for j_id, j_cont in jobs.items():
            if 'accession' in j_cont and j_cont['status'] != 'created':
                log=""
                if 'log' in j_cont:
                    log=j_cont['log']

                wdir = os.getcwd()
                log_path = wdir+"/"+log
                errorInfo = "Status: "
                #report log status
                if j_cont['status']== "job_finished":
                    errorInfo+="OK"
                    l.info(f"{j_cont['accession']}\t{j_cont['name']}\t{errorInfo}" )
                else:
                    errorInfo+="ERROR"
                    l.info(f"{j_cont['accession']}\t{j_cont['name']}\t{errorInfo}\t{log_path}" )

                if j_cont['status']== "job_error":
                    if j_cont['name']== "get_experiment_metadata":
                        if os.path.isfile(log_path):
                            file = open( log_path, "r")
                            err=[]
                            for line in file:                                                                                                                                                                            
                                if re.search("error", line, re.IGNORECASE):                                                                                                                                              
                                    err.append(line.strip())                                                                                                                                                               
                            if len(err)>0:                                                                                                                                                                                 
                                l.info(f"-- error/s found in get_experiment_metadata log: {len(err)}" )                                                                                                                           
                                for i in range(len(err)):                                                                                                                                                                  
                                    l.info(f"---{i+1}\t{j_cont['accession']}\t{err[i]}" )                                                                                                                                    
                                                                                                                                                                                                                         
                    if j_cont['name']== "produce_recalculations_call":                                                                                                                                                   
                        if os.path.isfile(log_path):                                                                                                                                                                     
                            file = open( log_path, "r")                                                                                                                                                                  
                            err0=[]                                                                                                                                                                                         
                            err1=[]                                                                                                                                                                                         
                            for line in file:                                                                                                                                                                            
                                if re.search("error in rule", line, re.IGNORECASE):                                                                                                                                                                                                                                                                                        
                                    err0.append(line.strip())                                                                                                                                                               
                                if re.search("for error message\)", line, re.IGNORECASE):                                                                                                                                                                                                                                                                                 
                                    err1.append(line.strip())                                                                                                                                                               
                            if len(err0)>0:                                                                                                                                                                                 
                                l.info(f"-- error in rules: {len(err0)}" )                                                                                                                                                  
                                l.info(f"-- error messages: {len(err1)}" )                                                                                                                                                  
                                seen=[]                                                                                                                                  
                                for i in range(len(err1)):
                                    seen.append(err1[i])                                                                                                                                                                  
                                    err=err1[i].replace("(check log file(s) for error message)","")                                                                                                                                                            
                                    l.info(f"---{i+1}\t{j_cont['accession']}\t{err0[i]}\t{err}\tattempt: {seen.count(err1[i])} ")                                                                                                                               
                                                                                                                                                                                                                                                      
                jobs_to_remove.append(j_id)
                #sys.stdout.flush()

        for j_id in jobs_to_remove:
            del jobs[j_id]

    #if msg.get('level') == 'progress':
    #    cur = msg['done']
    #    total = msg['total']
    #    print(get_bar(cur, total) + f" {cur}/{total}")
    #if 'msg' in msg and msg['msg'] == 'removed all locks':
    #    sys.stdout.close()
