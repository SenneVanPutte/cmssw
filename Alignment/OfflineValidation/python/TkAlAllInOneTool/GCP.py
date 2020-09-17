import os
import copy

def GCP(config, validationDir):
    ##List with all jobs
    jobs = []

    for comparison in config["validations"]["GCP"]["compare"]:
        for ali_pair in config["validations"]["GCP"]["compare"][comparison]:
            ref_name  = copy.deepcopy(config["validations"]["GCP"]["compare"][comparison][ali_pair]["reference"])
            comp_name = copy.deepcopy(config["validations"]["GCP"]["compare"][comparison][ali_pair]["compared"])
            IOVpair_list = []
            IOVali_list = []

            # Construct pairs from IOVlist
            IOV_list = []
            if "IOVlist" in config["validations"]["GCP"]["compare"][comparison][ali_pair]: IOV_list = copy.deepcopy(config["validations"]["GCP"]["compare"][comparison][ali_pair]["IOVlist"])
            IOV_list.sort()
            for idx,IOV in enumerate(IOV_list):
                if ref_name == comp_name:
                    IOV_pair = str(IOV)+'_vs_'+str(IOV_list[0])
                    IOV_ali_r = ref_name+'_'+str(IOV_list[0])
                    IOV_ali_c = comp_name+'_'+str(IOV)
                else:
                    IOV_pair = str(IOV)+'_vs_'+str(IOV)
                    IOV_ali_r = ref_name+'_'+str(IOV)
                    IOV_ali_c = comp_name+'_'+str(IOV)
                if IOV_pair not in IOVpair_list: IOVpair_list.append(IOV_pair)
                if IOV_ali_r not in IOVali_list: IOVali_list.append(IOV_ali_r)
                if IOV_ali_c not in IOVali_list: IOVali_list.append(IOV_ali_c)

            # Read explicit pairs from IOVpairs
            pair_list = []
            if "IOVpairs" in config["validations"]["GCP"]["compare"][comparison][ali_pair]: pair_list = copy.deepcopy(config["validations"]["GCP"]["compare"][comparison][ali_pair]["IOVpairs"])
            for IOV_p in pair_list:
                IOV_pair = str(IOV_p[0])+'_vs_'+str(IOV_p[1])
                IOV_ali_r = ref_name+'_'+str(IOV_p[1])
                IOV_ali_c = comp_name+'_'+str(IOV_p[0])
                if IOV_pair not in IOVpair_list: IOVpair_list.append(IOV_pair)
                if IOV_ali_r not in IOVali_list: IOVali_list.append(IOV_ali_r)
                if IOV_ali_c not in IOVali_list: IOVali_list.append(IOV_ali_c)

            # GCP Ntuple job preps
            for IOV_ali in IOVali_list:
                ali = IOV_ali.split('_')[0]
                IOV = int(IOV_ali.split('_')[1])
                workDir = "{}/GCP/{}/{}/{}".format(validationDir, comparison, 'Ntuples', IOV_ali)
                
                # local config 
                local = {}
                local["output"] = "{}/{}/{}/{}/{}".format(config["LFS"], config["name"], comparison, 'Ntuples', IOV_ali)
                local["aligment"] = copy.deepcopy(config["alignments"][ali])
                local["validation"] = {}
                local["validation"]['GCP'] = copy.deepcopy(config["validations"]["GCP"][comparison])
                local["validation"]['IOV'] = IOV

                # job info
                job = {
                    "name": "GCP_{}_Ntuple_{}".format(comparison, IOV_ali),
                    "dir": workDir,
                    "exe": "cmsRun",
                    "cms-config": "{}/src/Alignment/OfflineValidation/python/TkAlAllInOneTool/GCP_Ntuples_cfg.py".format(os.environ["CMSSW_BASE"]),
                    "run-mode": "Condor",
                    "dependencies": [],
                    "config": local, 
                }

                jobs.append(job)

            # Comparison job preps
            for IOV_pair in IOVpair_list:
                ref_IOV  = int(IOV_pair.split('_vs_')[1])
                comp_IOV = int(IOV_pair.split('_vs_')[0])
                # workdir for each GCP, alignment pair and IOV pair
                workDir = "{}/GCP/{}/{}/{}".format(validationDir, comparison, ali_pair, IOV_pair)
               
                # local config
                local = {} 
                local["output"] = "{}/{}/{}/{}/{}".format(config["LFS"], config["name"], comparison, ali_pair, IOV_pair)
                local["aligments"] = {}
                local["aligments"]["ref"]  = copy.deepcopy(config["alignments"][ref_name])
                local["aligments"]["comp"] = copy.deepcopy(config["alignments"][comp_name])
                local["validation"] = {}
                local["validation"]['GCP'] = copy.deepcopy(config["validations"]["GCP"][comparison])
                local["validation"]["IOVref"] = ref_IOV
                local["validation"]["ALIref"] = ref_name
                local["validation"]["IOVcomp"] = comp_IOV
                local["validation"]["ALIcomp"] = comp_name

                # job info
                job = {
                    "name": "GCP_{}_{}_{}".format(comparison, ali_pair, IOV_pair),
                    "dir": workDir,
                    "exe": "GCP",
                    "run-mode": "Condor",
                    "dependencies": [],
                    "config": local, 
                }

                # setup dependancies 
                for j in jobs:
                    if not comparison in j['name']: continue
                    if not 'Ntuple' in j['name']: continue
                    if ref_name in j['name'] and str(ref_IOV) in j['name']: 
                        job["dependencies"].append(j['name'])
                        local["input_ref"] = j['config']['output']
                    if comp_name in j['name'] and str(comp_IOV) in j['name']: 
                        job["dependencies"].append(j['name'])
                        local["input_comp"] = j['config']['output']

                jobs.append(job)

    return jobs
