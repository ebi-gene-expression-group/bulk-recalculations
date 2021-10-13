import sys
import logging
import os
import re

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

        #for k, v in msg.items():
        #    print(f"Key: {k}")
        #    print(f"Value: {v}")

        jobs_to_remove = []
        for j_id, j_cont in jobs.items():
            if 'accession' in j_cont and j_cont['status'] != 'created':
                log=""
                if 'log' in j_cont:
                    log=j_cont['log']

                wdir = os.getcwd()
                log_path = wdir+"/"+log
                errorInfo = "Status: "
                tail_e="(check log file(s) for error message)"

                if j_cont['status']== "job_finished":
                    errorInfo+="OK"
                else:
                    errorInfo+="ERROR"
                #report log
                l.info(f"{j_cont['accession']}\t{j_cont['name']}\t{j_cont['status']}\t{log_path}\t {errorInfo}" )

                if j_cont['status']== "job_error":
                    if j_cont['name']== "get_experiment_metadata":
                        if os.path.isfile(log_path):
                            file = open( log_path, "r")
                            err=[]
                            for line in file:                                                                                                                                                                            
                                if re.search("error", line, re.IGNORECASE):                                                                                                                                              
                                    err.append(line.strip())                                                                                                                                                               
                            if len(err)>0:                                                                                                                                                                                 
                                l.info(f"-- error/s found get_experiment_metadata: {len(err)}" )                                                                                                                           
                                for i in range(len(err)):                                                                                                                                                                  
                                    l.info(f"---{i}\t{j_cont['accession']}\t{err[i]}" )                                                                                                                                    
                                                                                                                                                                                                                         
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
                                l.info(f"-- restart-times: {(len(err1)/len(set(err1)))-1 } ")                                                                                                                                  
                                for i in range(len(err0)):                                                                                                                                                                  
                                    err=err1[i].replace(tail_e,"")                                                                                                                                                            
                                    l.info(f"---{i}\t{j_cont['accession']}\t{err0[i]}\t{err} ")                                                                                                                               
                                                                                                                                                                                                                                                      
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
