import yaml
import csv
import math
import ast

Carbon = 12.00000
Hydrogen = 1.0078250350
Nitrogen = 14.0030740000
Oxygen = 15.9949146300
Sulfur = 31.9720707000
Phosphorus = 30.9737620000
Iodine = 126.904473
Chlorine = 34.9688527300
Fluorine = 18.9984032200
Bromine = 78.918361000
Sodium = 22.989767700
Potassium = 38.9637069000
Calcium = 39.9625912000
Electron = 0.000549
Water = 18.0105647000
Proton = 1.007276035

ms1_tolerance = 0.00001
ms2_tolerance = 0.00002

def match_frags():
    
    print ("input gly1 results with theoretical frags")

    filename = raw_input()
    f_csv = csv.DictReader(open(filename))

    print("enter data filename")
    datafile = raw_input()    
    stream = open(datafile,'r')
    data = []
    data.append(yaml.load(stream))

    results = []

    for lines in f_csv: # reading each line in theoretical ions file generated from gly1 results
        
        for rows in data: #reading yaml data scan by scan
                
            for more in rows['peaks']:

                real_ions = []
                real_ions = more['mass']  #all ions in scan assigned to list real_ions

                #print ("intensity array", more['intensity'])
                #print ("number of elements" , more['num'])

                info = more['scans'][:1]  #this takes the information from the first scan in case there are multiple merged scans
                for num in info:
                    
                    #print ("mz", info['mz'])
                    charge = num['z']
                    mass = num['mz']
                    ntr_mass = ((mass * charge)-(charge * Proton)) 
                    #print ("ntr-mass", ntr_mass)
                    
                    obs_mass = float(lines['Obs_Mass'])
                    ppm = ((ntr_mass-obs_mass)/ntr_mass)
                
                if math.fabs(ppm) <= ms1_tolerance:
                    Oxonium_ions = ast.literal_eval(lines['Oxonium_ions'])

                    oxoniums =[]
                    for ox_ion in Oxonium_ions:
                        for ob_ions in real_ions:                        
                            oxonium_ppm = (((ob_ions+Proton) - ox_ion)/(ob_ions+Proton))
                            if math.fabs(oxonium_ppm) <= ms2_tolerance:
                                oxoniums.append({"ion":(ob_ions+Proton), "ppm_error":oxonium_ppm, "key":ox_ion})
                                #print ("\n", ob_ions, "ppm error:", oxonium_ppm) 
                            else:
                                pass
                    oxo_found = len(oxoniums)

                    #checking for b and y ions and ions with HexNAc:
                    
                    b_ions = ast.literal_eval(lines['bare_b_ions'])
                    b_len = float(len(b_ions))
                    #print "\n, b ions:"
                    b_type = []
                    all_b_ions = []
                    for theo_ions in b_ions:
                        frag_mass = float(theo_ions["mass"])
                        prot_ion = frag_mass-Hydrogen #change this to proton once the error is found in Han's code
                        for obs_ions in real_ions:
                            tandem_ppm = float((obs_ions-(prot_ion))/obs_ions)
                            if math.fabs(tandem_ppm) <= ms2_tolerance:
                                b_type.append({"obs_ion":(obs_ions+Proton), "key":theo_ions["key"], "ppm_error":tandem_ppm})
                                all_b_ions.append({"obs_ion":(obs_ions+Proton), "key":theo_ions["key"], "ppm_error":tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm)
                            else:
                                pass

                    len_found = float(len(b_type))
                    if len_found != 0: 
                        b_found = (len_found/b_len)*100
                    else:
                        b_found = 0

                        
                    b_HexNAc_ions = ast.literal_eval(lines['b_ions_with_HexNAc'])
                    b_HexNAc_len = float(len(b_HexNAc_ions))
                    #print "\n, b ions:"
                    b_HexNAc_type = []
                    for theo_ions in b_HexNAc_ions:
                        frag_mass = float(theo_ions["mass"])
                        prot_ion = frag_mass-Hydrogen #change this to proton once the error is found in Han's code
                        for obs_ions in real_ions:
                            tandem_ppm = float((obs_ions-(prot_ion))/obs_ions)
                            if math.fabs(tandem_ppm) <= ms2_tolerance:
                                b_HexNAc_type.append({"obs_ion":(obs_ions+Proton), "key":theo_ions["key"], "ppm_error":tandem_ppm})
                                all_b_ions.append({"obs_ion":(obs_ions+Proton), "key":theo_ions["key"].split("+")[0], "ppm_error":tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm)
                            else:
                                pass
                    len_found = float(len(b_HexNAc_type))
                    if len_found != 0: 
                        b_HexNAc_found = (len_found/b_HexNAc_len)*100
                    else:
                        b_HexNAc_found = 0
                    



                    y_ions = list(ast.literal_eval(lines['bare_y_ions']))
                    #print "\n, y ions:"
                    y_len = float(len(y_ions))
                    y_type= []
                    all_y_ions = []
                    for theo_ions in y_ions:
                        frag_mass = float(theo_ions["mass"])
                        prot_ion = frag_mass-Hydrogen #change this to proton once the error is found in Han's code
                        for obs_ions in real_ions:
                            tandem_ppm = float((obs_ions-(prot_ion))/obs_ions)
                            if math.fabs(tandem_ppm) <= ms2_tolerance:
                                y_type.append({"obs_ion":(obs_ions+Proton), "key":theo_ions["key"], "ppm_error":tandem_ppm})
                                all_y_ions.append({"obs_ion":(obs_ions+Proton), "key":theo_ions["key"], "ppm_error":tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm, "\n")
                            else:
                                pass
                    len_found = float(len(y_type))
                    if len_found != 0:
                        y_found = (len_found/y_len)*100
                    elif len_found ==0:
                        y_found = 0
                  

                    y_HexNAc_ions = ast.literal_eval(lines['y_ions_with_HexNAc'])
                    y_HexNAc_len = float(len(y_HexNAc_ions))
                    #print "\n, b ions:"
                    y_HexNAc_type = []
                    for theo_ions in y_HexNAc_ions:
                        frag_mass = float(theo_ions["mass"])
                        prot_ion = frag_mass-Hydrogen #change this to proton once the error is found in Han's code
                        for obs_ions in real_ions:
                            tandem_ppm = float((obs_ions-(prot_ion))/obs_ions)
                            if math.fabs(tandem_ppm) <= ms2_tolerance:
                                y_HexNAc_type.append({"obs_ion":(obs_ions+Proton), "key":theo_ions["key"], "ppm_error":tandem_ppm})
                                all_y_ions.append({"obs_ion":(obs_ions+Proton), "key":theo_ions["key"].split("+")[0], "ppm_error":tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm)
                            else:
                                pass
                    len_found = float(len(y_HexNAc_type))
                    if len_found != 0: 
                        y_HexNAc_found = (len_found/y_HexNAc_len)*100
                    else:
                        y_HexNAc_found = 0

                        
                    #checking for stub ions

                    stub_ions = ast.literal_eval(lines['pep_stub_ions'])
                    stub_type = []
                    for theo_ions in stub_ions:
                        frag_mass = float(theo_ions["mass"])
                        prot_ion = frag_mass-Proton
                        for obs_ions in real_ions:
                            tandem_ppm = float((obs_ions-(prot_ion))/obs_ions)
                            if math.fabs(tandem_ppm) <= ms2_tolerance:
                                stub_type.append({"obs_ion":(obs_ions+Proton), "key":theo_ions["key"], "ppm_error":tandem_ppm})
                                #print (obs_ions, theo_ions["mass"], theo_ions["key"], tandem_ppm, "\n")
                            else:
                                pass
                    if len(stub_type) != 0:
                        stub_found = len (stub_type)
                    else:
                        stub_found = 0


                    results.append({"MS1_Score":lines["MS1_Score"],"Obs_Mass":lines["Obs_Mass"],"Calc_mass":lines["Calc_mass"],"ppm_error":lines["ppm_error"],"Peptide":lines["Peptide"],"Peptide_mod":lines["Peptide_mod"],"Glycan":lines["Glycan"],"vol":lines["vol"],"glyco_sites":lines["glyco_sites"],"startAA":lines["startAA"],"endAA":lines["endAA"],"Seq_with_mod":lines["Seq_with_mod"],"Glycopeptide_identifier":lines["Seq_with_mod"] + lines["Glycan"],"Oxonium_ions":oxoniums,"#_Oxonium_found":oxo_found,"bare_b_ions":b_type,"%_b_ions_found":b_found,"bare_y_ions":y_type,"%_y_ions_found":y_found, "b_ions_with_HexNAc":b_HexNAc_type,"%_b_ions_w/_HexNAc_found":b_HexNAc_found,"y_ions_with_HexNAC":y_HexNAc_type,"%_y_ions_w/_HexNAc_found":y_HexNAc_found,"b_ion_coverage":all_b_ions,"y_ion_coverage":all_y_ions,"Stub_ions":stub_type,"#_stubs_found":stub_found})
                
                        
                else:
                    pass

    
    keys = ["MS1_Score","Obs_Mass","Calc_mass","ppm_error","Peptide","Peptide_mod","Glycan","vol","glyco_sites","startAA","endAA","Seq_with_mod","Glycopeptide_identifier","Oxonium_ions","#_Oxonium_found","bare_b_ions","%_b_ions_found","bare_y_ions","%_y_ions_found","b_ions_with_HexNAc","%_b_ions_w/_HexNAc_found","y_ions_with_HexNAC","%_y_ions_w/_HexNAc_found","b_ion_coverage","y_ion_coverage","Stub_ions","#_stubs_found"]
    print("enter output filename")
    filename = raw_input()
    f = open(filename,'wb')
    dict_writer = csv.DictWriter(f,keys)
    dict_writer.writer.writerow(keys)
    dict_writer.writerows(results)                     
match_frags()