import pymzml
import matplotlib.pyplot as plt
import os
import sqlite3
import struct
import csv


file_path = 'names_identifiers.tsv'
data = {}

with open(file_path, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        print(row)
        try:
            key, value = row[0], row[1]
            data[key] = value
        except IndexError:
            print(row[0])




conn = sqlite3.connect("msms_db.db")
conn.execute('CREATE TABLE IF NOT EXISTS MSMS (id INTEGER PRIMARY KEY, name TEXT, identifier TEXT, adduct TEXT, cid TEXT, window TEXT, tof TEXT, mz TEXT, spectra BLOB)')


def get_mzml_files_in_subfolders(folder_path):
    mzml_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".mzML"):
                mzml_files.append(os.path.join(root, file))
    return mzml_files


def list_to_blob(input_list):
    return b''.join(struct.pack('dq', float_val, int(int_val)) for float_val, int_val in input_list)

def insert_item_in_sqlite(mzml_file, input_list):
    if "FMP10" in mzml_file.split("\\")[-2] or "FMP10" in mzml_file.split("/")[-2]:
        name = mzml_file.split("\\")[-2]
        identifier = data.get(name, "")
        filename = mzml_file.split("\\")[-1]
        print(mzml_file)
        
        window = filename.split("_")[3]
        cid = filename.split("_")[6]
        tof = ".".join(filename.split("_")[7:9]).strip("tof")
        mz = ".".join(filename.split("_")[1:3]).strip("mz")
        adduct = mz
        blob_data = list_to_blob(input_list)
        print("name:", name, "adduct:", adduct, "mz:", mz, "cid:", cid, "tof:", tof, "window:", window)
        conn.execute('INSERT INTO MSMS (name, identifier, adduct, cid, window, tof, mz, spectra) VALUES (?, ?, ?, ?, ?, ?, ?, ?)', (name, identifier, adduct, cid, window, tof, mz, blob_data))
        conn.commit()
    else:

        try:
            name = mzml_file.split("\\")[-2]
            identifier = data.get(name, None)
            filename = mzml_file.split("\\")[-1]
            print(mzml_file)
            adduct = filename.split("_")[3]
            window = filename.split("_")[4]
            cid = filename.split("_")[7]
            tof = ".".join(filename.split("_")[8:10]).strip("tof")
            mz = ".".join(filename.split("_")[1:3]).strip("mz")
            blob_data = list_to_blob(input_list)
            #print("name:", name, "adduct:", adduct, "mz:", mz)
            conn.execute('INSERT INTO MSMS (name, identifier, adduct, cid, window, tof, mz, spectra) VALUES (?, ?, ?, ?, ?, ?, ?, ?)', (name, identifier, adduct, cid, window, tof, mz, blob_data))
            conn.commit()
        
        except IndexError:
            name = mzml_file.split("/")[-2]
            filename = mzml_file.split("/")[-1]
            print(mzml_file)
            adduct = filename.split("_")[3]
            window = filename.split("_")[4]
            cid = filename.split("_")[7]
            tof = ".".join(filename.split("_")[8:10]).strip("tof")
            mz = ".".join(filename.split("_")[1:3]).strip("mz")
            blob_data = list_to_blob(input_list)
            #print("name:", name, "adduct:", adduct, "mz:", mz)
            conn.execute('INSERT INTO MSMS (name, identifier, adduct, cid, window, tof, mz, spectra) VALUES (?, ?, ?, ?, ?, ?, ?, ?)', (name, identifier, adduct, cid, window, tof, mz, blob_data))
            conn.commit()
        except: 
            print(name)
    
    

folder_path = '//mms2.farmbio.uu.se/imaging3/NewProjects2023/Metabolite_Identification/202306_FMP10_MSMS_standards'
mzml_files = get_mzml_files_in_subfolders(folder_path)

print(f"Found {len(mzml_files)} .mzML files:")

for mzml_file in mzml_files:
    run = pymzml.run.Reader(mzml_file)
    for n, spec in enumerate(run):
        spectra = spec.remove_noise(noise_level=spec.estimated_noise_level()*1.5).peaks("centroided")       
        insert_item_in_sqlite(mzml_file, spectra)
            

        
        




















