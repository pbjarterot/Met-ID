import { save, open } from '@tauri-apps/api/dialog';
import { invoke } from '@tauri-apps/api';

async function read_input_csv(result: string, delimiter: string, transpose: boolean, deleted_rows: number[]) {
  const parsed_csv: string[][] = await invoke("read_mass_error_csv", {path: result, delimiter: delimiter, transpose: transpose, deletedrows: deleted_rows});
  return parsed_csv
}

export async function convertTableToCSV(table_id: string) {
  const table = document.getElementById(table_id) as HTMLTableElement;
  const rows = table.querySelectorAll('tr');

  let csvContent: string;
  if (table_id === "ms1-datatable") {
    csvContent = ms1_table_content(rows);
    
  } else {
    csvContent = mass_error_table_content(rows);
  }
  try {
    // Open a save file dialog
    const result = await save({
      filters: [{name: 'CSV Files', extensions: ["csv",]}],
      defaultPath: 'met-id-results.csv', // Default file name
    });

    await invoke("save_csv", {path:result, csvcontent: csvContent});

  } catch (error) {
    console.error('Error saving CSV file:', error);
  }
}

function mass_error_table_content(rows:NodeListOf<HTMLTableRowElement>) {
  let csvContent = '';
  rows.forEach((row) => {
    const rowData: string[] = [];
    const cells = row.querySelectorAll('td, th');

    // Loop through only the first two cells (columns)
    for (let i = 0; i < Math.min(cells.length, 2); i++) {
      rowData.push(`${String(cells[i].textContent).trim()}`);
    }

    csvContent += rowData.join(',') + '\n';
  });
  return csvContent
}

function countBrTagsInCell(cell: HTMLTableCellElement): number {
  // Use querySelectorAll to select all <br> tags within the cell
  const brTags = cell.querySelectorAll('br');
  // Return the count of <br> tags
  return brTags.length;
}

function ms1_table_content(rows:NodeListOf<HTMLTableRowElement>) {
  let csvContent = '';
  rows.forEach((row) => {
    
    const cells = row.querySelectorAll('td, th');
    if (cells[3]) {
      if (cells.length > 0 && cells[3].textContent != "" ) {
        const brCount = countBrTagsInCell(cells[3] as HTMLTableCellElement);
        for (let br_idx = 0; br_idx < brCount + 1; br_idx++) {
			const rowData: string[] = [];
            for (let i = 1; i < cells.length-2; i++) {
                const cellContent = cells[i].innerHTML.split(/<br>/g)[br_idx];
                rowData.push(`${String(cellContent)}`);
            }
            csvContent += rowData + "\n";
        }
    }
    }
    
  });
  return csvContent
}

export async function importFromCSV() {
  let result: string = await open({directory: false, multiple: false, filters: [{name: 'Delimiter Separated Files', extensions:['csv', 'tsv']}, {name: 'Text Files', extensions:['txt']}]}) as string;
  let file_contents: string[][] = await read_input_csv(result as string, ",", false, [0, 1] );
  
  return file_contents
}


function zip<T>(...arrays: T[][]): T[][] {
  const length = Math.min(...arrays.map((arr) => arr.length));
  return Array.from({ length }, (_, i) => arrays.map((arr) => arr[i]));
}

export async function saveCsv(oMass: string[], aMass: string[]): Promise<void> {
  try {
    oMass.unshift("Observed Mass");
    aMass.unshift("Adjusted Mass");

    const csvContent = zip(oMass, aMass).join('\n');

    // Open a save file dialog
    const result = await save({
      filters: [{name: 'CSV Files', extensions: [".csv",]}],
      defaultPath: 'adjusted_masses.csv', // Default file name
    });

    await invoke("save_csv", {path:result, csvcontent: csvContent});


  } catch (error) {
    console.error('Error saving CSV file:', error);
  }
}