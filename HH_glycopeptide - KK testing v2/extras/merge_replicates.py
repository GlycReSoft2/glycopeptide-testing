import csv
import math

def merge_output():
    print "enter scored csv output filename #1"
    csv_1= raw_input()
    temp_1 = csv.DictReader(open(csv_1))
    f_csv_1 = []
    for row in temp_1:
        f_csv_1.append(row)
    print "enter scored csv output filename #2"
    csv_2= raw_input()
    temp_2 = csv.DictReader(open(csv_2))
    f_csv_2 = []
    for row in temp_2:
        f_csv_2.append(row)
    print "enter scored csv output filename #3"
    csv_3= raw_input()
    temp_3 = csv.DictReader(open(csv_3))
    f_csv_3 = []
    for row in temp_3:
        f_csv_3.append(row)
    new_csv = [] 

    
    for i in f_csv_1:
        
        pep_id = i["Glycopeptide_identifier"]
        
        for j in f_csv_2:
            if i["Glycopeptide_identifier"] == j["Glycopeptide_identifier"]:
                for k in f_csv_3:
                   if j["Glycopeptide_identifier"] == k["Glycopeptide_identifier"]:
                                            
                       average_MS1_Score = float ((float(i["MS1_Score"])+float(j["MS1_Score"])+float(k["MS1_Score"]))/3)
                       MS1_Score_std_dev = math.sqrt((((float(i["MS1_Score"])-average_MS1_Score)**2)+((float(j["MS1_Score"])-average_MS1_Score)**2)+((float(k["MS1_Score"])-average_MS1_Score)**2))/3)

                       avg_b_HexNAc_coverage = float ((float(i["%b_ion_with_HexNAc_coverage"])+float(j["%b_ion_with_HexNAc_coverage"])+float(k["%b_ion_with_HexNAc_coverage"]))/3)
                       b_HexNAc_coverage_std_dev = math.sqrt((((float(i["%b_ion_with_HexNAc_coverage"])-avg_b_HexNAc_coverage)**2)+((float(j["%b_ion_with_HexNAc_coverage"])-avg_b_HexNAc_coverage)**2)+((float(k["%b_ion_with_HexNAc_coverage"])-avg_b_HexNAc_coverage)**2))/3)

                       avg_y_HexNAc_coverage = float ((float(i["%y_ion_with_HexNAc_coverage"])+float(j["%y_ion_with_HexNAc_coverage"])+float(k["%y_ion_with_HexNAc_coverage"]))/3)
                       y_HexNAc_coverage_std_dev = math.sqrt((((float(i["%y_ion_with_HexNAc_coverage"])-avg_y_HexNAc_coverage)**2)+((float(j["%y_ion_with_HexNAc_coverage"])-avg_y_HexNAc_coverage)**2)+((float(k["%y_ion_with_HexNAc_coverage"])-avg_y_HexNAc_coverage)**2))/3)

                       avg_b_coverage = float ((float(i["%percent_b-ion_coverage"])+float(j["%percent_b-ion_coverage"])+float(k["%percent_b-ion_coverage"]))/3)
                       b_coverage_std_dev = math.sqrt((((float(i["%percent_b-ion_coverage"])-avg_b_coverage)**2)+((float(j["%percent_b-ion_coverage"])-avg_b_coverage)**2)+((float(k["%percent_b-ion_coverage"])-avg_b_coverage)**2))/3)

                       avg_y_coverage = float ((float(i["%percent_y-ion_coverage"])+float(j["%percent_y-ion_coverage"])+float(k["%percent_y-ion_coverage"]))/3)
                       y_coverage_std_dev = math.sqrt((((float(i["%percent_y-ion_coverage"])-avg_y_coverage)**2)+((float(j["%percent_y-ion_coverage"])-avg_y_coverage)**2)+((float(k["%percent_y-ion_coverage"])-avg_y_coverage)**2))/3)

                       avg_stubs = float((float(i["#_of_stubs_found"])+float(j["#_of_stubs_found"])+float(k["#_of_stubs_found"]))/3)
                       std_dev_stubs = math.sqrt((((float(i["#_of_stubs_found"])-avg_stubs)**2)+((float(j["#_of_stubs_found"])-avg_stubs)**2)+((float(k["#_of_stubs_found"])-avg_stubs)**2))/3)

                       avg_MS2_score = float ((float(i["MS2_score"])+float(j["MS2_score"])+float(k["MS2_score"]))/3)
                       MS2_score_std_dev = math.sqrt((((float(i["MS2_score"])-avg_MS2_score)**2)+((float(j["MS2_score"])-avg_MS2_score)**2)+((float(k["MS2_score"])-avg_MS2_score)**2))/3)
                       avg_vol = float ((float(i["vol"])+float(j["vol"])+float(k["vol"]))/3)
                       vol_std_dev = math.sqrt((((float(i["vol"])-avg_vol)**2)+((float(j["vol"])-avg_vol)**2)+((float(k["vol"])-avg_vol)**2))/3)

                       new_csv.append({"average_MS1_Score":average_MS1_Score, "std_dev_MS1_score":MS1_Score_std_dev,"Obs_Mass":i["Obs_Mass"],"Calc_mass":i["Calc_mass"],"ppm_error":i["ppm_error"],"Peptide":i["Peptide"],"Peptide_mod":i["Peptide_mod"],"Glycan":i["Glycan"],"avg_vol":avg_vol,"vol_std_dev":vol_std_dev,"glyco_sites":i["glyco_sites"],"startAA":i["startAA"],"endAA":i["endAA"],"Seq_with_mod":i["Seq_with_mod"],"Glycopeptide_identifier":i["Glycopeptide_identifier"],"avg_%percent_b-ion_with_HexNAc_coverage":avg_b_HexNAc_coverage,"std_dev_b_ions_with_HexNAc_coverage":b_HexNAc_coverage_std_dev,"avg_%percent_y-ion_with_HexNAc_coverage":avg_y_HexNAc_coverage,"std_dev_y_ions_with_HexNAc_coverage":y_HexNAc_coverage_std_dev,"avg_%percent_b-ion_coverage":avg_b_coverage,"std_dev_b_coverage":b_coverage_std_dev, "avg_%percent_y-ion_coverage":avg_y_coverage, "std_dev_y_coverage": y_coverage_std_dev,"avg_#_of_stubs_found":avg_stubs, "std_dev_subs":std_dev_stubs ,"avg_MS2_score":avg_MS2_score, "std_dev_MS2_score":MS2_score_std_dev})
                
                   else:
                       pass
            else:
                pass
                        
    keys = ["average_MS1_Score", "std_dev_MS1_score","Obs_Mass","Calc_mass","ppm_error","Peptide","Peptide_mod","Glycan","avg_vol","vol_std_dev","glyco_sites","startAA","endAA","Seq_with_mod","Glycopeptide_identifier", "avg_%percent_b-ion_with_HexNAc_coverage","std_dev_b_ions_with_HexNAc_coverage","avg_%percent_y-ion_with_HexNAc_coverage","std_dev_y_ions_with_HexNAc_coverage","avg_%percent_b-ion_coverage","std_dev_b_coverage","avg_%percent_y-ion_coverage", "std_dev_y_coverage","avg_#_of_stubs_found", "std_dev_subs" ,"avg_MS2_score", "std_dev_MS2_score"]

    print("enter output filename")
    filename = raw_input()
    f = open(filename,'wb')
    dict_writer = csv.DictWriter(f,keys)
    dict_writer.writer.writerow(keys)
    dict_writer.writerows(new_csv)
    f.close()
    
merge_output()