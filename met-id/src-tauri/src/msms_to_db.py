
import sqlite3
import sys
import struct
import pymzml
import urllib.request


def list_to_blob(input_list):
    return b''.join(struct.pack('dq', float_val, int(int_val)) for float_val, int_val in input_list)


def insert_item_in_sqlite(input_list, name, adduct, mz, cid, tof, window, identifier, database_path):
	conn = sqlite3.connect(database_path)
	conn.execute('CREATE TABLE IF NOT EXISTS MSMS (id INTEGER PRIMARY KEY, name TEXT, identifier TEXT, adduct TEXT, cid TEXT, window TEXT, tof TEXT, mz TEXT, spectra BLOB)')

	blob_data = list_to_blob(input_list)
	print("name:", name, "adduct:", adduct, "mz:", mz, "cid:", cid, "tof:", tof, "window:", window)
	conn.execute('INSERT INTO MSMS (name, identifier, adduct, cid, window, tof, mz, spectra) VALUES (?, ?, ?, ?, ?, ?, ?, ?)', (name, identifier, adduct, cid, window, tof, mz, blob_data))
	conn.commit()


if __name__ == "__main__":
	print("here1")
	mzml_file = sys.argv[1]
	name = sys.argv[2]
	adduct = sys.argv[3]
	mz = sys.argv[4]
	cid = sys.argv[5]
	tof = sys.argv[6]
	window = sys.argv[7]
	identifier = sys.argv[8]
	database_path = sys.argv[9]

	print("here")
	run = pymzml.run.Reader(mzml_file)
	for n, spec in enumerate(run):
		print(n, spec)
		spectra = spec.remove_noise(noise_level=spec.estimated_noise_level()*1.5).peaks("centroided")       
		insert_item_in_sqlite(spectra, name, adduct, mz, cid, tof, window, identifier, database_path)