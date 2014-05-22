#! /usr/bin/python
import math

def MergeRows(SourceData):
	SourceData = sorted(SourceData, key=lambda k: k['Glycopeptide_identifier'])
	Storage = SourceData[0]
	#MergedList is used to store our final result, the merged list.
	MergedList = []
	for row in range(len(SourceData)):
		#newStorage is used to store the rows that have the same Glycopeptide identifier. Redeclaring it here to remove data from the last Glycopep ID.
		newStorage = []
		#input the first row.
		newStorage.append(Storage)
		if SourceData[row]['Glycopeptide_identifier'] == newStorage[0]['Glycopeptide_identifier']:
			newStorage.append(SourceData[row])
			#If this is the last row in SourceData, merge the rows now, instead of doing it in the "else" statement.
			if (row + 1) == len(SourceData):
				#Declaring variables to store N P R T V X Y Z columns.
				oi = []
				bbi = []
				byi = []
				biwh = []
				yiwh = []
				bic = []
				yic = []
				si = []
				#newrow is used to store our combined row
				newrow = []
				#input initial variables into newrow.
				newrow = newStorage[0]
				for row2 in newStorage:
					oi.append(row2["Oxonium_ions"])
					bbi.append(row2["bare_b_ions"])
					byi.append(row2["bare_y_ions"])
					biwh.append(row2["b_ions_with_HexNAc"])
					yiwh.append(row2["y_ions_with_HexNAC"])	
					bic.append(row2["b_ion_coverage"])
					yic.append(row2["y_ion_coverage"])
					si.append(row2["Stub_ions"])
				#Change the N P R T V X Y Z columns' values in newrow.
				newrow["Oxonium_ions"] = Mergedicts(oi)
				newrow["bare_b_ions"] = Mergedicts(bbi)
				newrow["bare_y_ions"] = Mergedicts(byi)
				newrow["b_ions_with_HexNAc"] = Mergedicts(biwh)
				newrow["y_ions_with_HexNAC"] = Mergedicts(yiwh)
				newrow["b_ion_coverage"] = Mergedicts(bic)
				newrow["y_ion_coverage"] = Mergedicts(yic)
				newrow["Stub_ions"] = Mergedicts(si)
				#Now put the newrow into the final output
				MergedList.append(newrow)
			
				#Declaring variables to store N P R T V X Y Z columns.
				oi = []
				bbi = []
				byi = []
				biwh = []
				yiwh = []
				bic = []
				yic = []
				si = []
				#newrow is used to store our combined row
				newrow = []
				#input initial variables into newrow.
				newrow = newStorage[0]
				for row2 in newStorage:
					oi.append(row2["Oxonium_ions"])
					bbi.append(row2["bare_b_ions"])
					byi.append(row2["bare_y_ions"])
					biwh.append(row2["b_ions_with_HexNAc"])
					yiwh.append(row2["y_ions_with_HexNAC"])	
					bic.append(row2["b_ion_coverage"])
					yic.append(row2["y_ion_coverage"])
					si.append(row2["Stub_ions"])
				#Change the N P R T V X Y Z columns' values in newrow.
				newrow["Oxonium_ions"] = Mergedicts(oi)
				newrow["bare_b_ions"] = Mergedicts(bbi)
				newrow["bare_y_ions"] = Mergedicts(byi)
				newrow["b_ions_with_HexNAc"] = Mergedicts(biwh)
				newrow["y_ions_with_HexNAC"] = Mergedicts(yiwh)
				newrow["b_ion_coverage"] = Mergedicts(bic)
				newrow["y_ion_coverage"] = Mergedicts(yic)
				newrow["Stub_ions"] = Mergedicts(si)
				#Now put the newrow into the final output
				MergedList.append(newrow)
		else:
			#The new row has another glycopep identifier, so we are storing it in Storage for the next loop. 
			Storage = SourceData[row]
			#Lets start merging the rows now.
			#Declaring variables to store N P R T V X Y Z columns.
			oi = []
			bbi = []
			byi = []
			biwh = []
			yiwh = []
			bic = []
			yic = []
			si = []
			#newrow is used to store our combined row
			newrow = []
			#input initial variables into newrow.
			newrow = newStorage[0]
			for row2 in newStorage:
				oi.append(row2["Oxonium_ions"])
				bbi.append(row2["bare_b_ions"])
				byi.append(row2["bare_y_ions"])
				biwh.append(row2["b_ions_with_HexNAc"])
				yiwh.append(row2["y_ions_with_HexNAC"])	
				bic.append(row2["b_ion_coverage"])
				yic.append(row2["y_ion_coverage"])
				si.append(row2["Stub_ions"])
			#Change the N P R T V X Y Z columns' values in newrow.
			newrow["Oxonium_ions"] = Mergedicts(oi)
			newrow["bare_b_ions"] = Mergedicts(bbi)
			newrow["bare_y_ions"] = Mergedicts(byi)
			newrow["b_ions_with_HexNAc"] = Mergedicts(biwh)
			newrow["y_ions_with_HexNAC"] = Mergedicts(yiwh)
			newrow["b_ion_coverage"] = Mergedicts(bic)
			newrow["y_ion_coverage"] = Mergedicts(yic)
			newrow["Stub_ions"] = Mergedicts(si)
			#Now put the newrow into the final output
			MergedList.append(newrow)
		
			#Declaring variables to store N P R T V X Y Z columns.
			oi = []
			bbi = []
			byi = []
			biwh = []
			yiwh = []
			bic = []
			yic = []
			si = []
			#newrow is used to store our combined row
			newrow = []
			#input initial variables into newrow.
			newrow = newStorage[0]
			for row2 in newStorage:
				oi.append(row2["Oxonium_ions"])
				bbi.append(row2["bare_b_ions"])
				byi.append(row2["bare_y_ions"])
				biwh.append(row2["b_ions_with_HexNAc"])
				yiwh.append(row2["y_ions_with_HexNAC"])	
				bic.append(row2["b_ion_coverage"])
				yic.append(row2["y_ion_coverage"])
				si.append(row2["Stub_ions"])
			#Change the N P R T V X Y Z columns' values in newrow.
			newrow["Oxonium_ions"] = Mergedicts(oi)
			newrow["bare_b_ions"] = Mergedicts(bbi)
			newrow["bare_y_ions"] = Mergedicts(byi)
			newrow["b_ions_with_HexNAc"] = Mergedicts(biwh)
			newrow["y_ions_with_HexNAC"] = Mergedicts(yiwh)
			newrow["b_ion_coverage"] = Mergedicts(bic)
			newrow["y_ion_coverage"] = Mergedicts(yic)
			newrow["Stub_ions"] = Mergedicts(si)
			#Now put the newrow into the final output
			MergedList.append(newrow)
			
	return MergedList
			
def Mergedicts(list):
	#"firstrow" stores the data we will output. Items from other rows are incorporated into it if they have a ppm error closer to 0, or if they have a row not existing in the first row.
	#in case the first list is zero, we need a loop here. "allzero" is used to check if all "rows in list" are empty.
	allzero = True
	for row in list:
		if len(row) != 0:
			allzero = False
			firstrow = row
			break;
	
	if allzero:
	#if the list is empty anyway, just return the same list back.
		return list
	else:
		#if the list isn't empty, merge them.		
		#Each of the "row in list" should be something like this: [{'ppm_error': -3.3899212497496274e-06, 'key': 'B3', 'obs_ion': 372.150116035}, {'ppm_error': -1.2204514502310785e-05, 'key': 'B4', 'obs_ion': 532.175476035}, {'ppm_error': 1.812995934605501e-05, 'key': 'B4', 'obs_ion': 532.191589035}, {'ppm_error': 1.4021538362178224e-07, 'key': 'B5', 'obs_ion': 645.266113035}, {'ppm_error': 4.080114450361758e-06, 'key': 'B7', 'obs_ion': 922.376038035}, {'ppm_error': 1.0112467130843996e-05, 'key': 'B14', 'obs_ion': 1741.8025510349999}]
		for row in list:
		#each of the "ReadingItems and firstrowReadingItems" should be something like this: {'ppm_error': -3.3899212497496274e-06, 'key': 'B3', 'obs_ion'}: 372.150116035}. If not, plz tell me.
			for ReadingItems in row:
				same_key_exists = False
				for firstrowReadingItems in range(len(firstrow)):
					if ReadingItems['key'] == firstrow[firstrowReadingItems]['key']:
						#if the same key exists, merge them.
						same_key_exists = True
						if math.fabs(float(ReadingItems['ppm_error'])) < math.fabs(float(firstrow[firstrowReadingItems][ppm_error])):
							firstrow[firstrowReadingItems] = ReadingItems
							break
				#if the row doesn't exist in our final answer, add it.
				if same_key_exists:
					same_key_exists = False
				else:
					firstrow.append(ReadingItems)
					
	return firstrow
