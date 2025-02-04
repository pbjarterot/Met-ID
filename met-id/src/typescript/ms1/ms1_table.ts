window.addEventListener("DOMContentLoaded", () => {
    document.getElementById("append-to-ms1-table")?.addEventListener("click", () => appendToMSTable());
    document.getElementById("ms1-table-delete-row")?.addEventListener("click", () => deleteRows())
    document.getElementById("ms1-table-head-check")?.addEventListener("click", () => selectRows())

    document.getElementById("toggle-show-identified")?.addEventListener("click", () => toggleRows())
});

let table_body = document.getElementById("ms1-table-body");

export function createTableRowHTML(input_mz: number): string {
  return `
      <tr>
          <td><input type="checkbox" /></td>
          <td class="clickable-cell">${input_mz}</td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"></td>
          <td class="clickable-cell"><span class="unavailable"></span></td>
      </tr>`;
}

function appendToMSTable() {
    let input_mass_element = document.getElementById("add-to-ms1-table-input") as HTMLInputElement;
    let input_mass = parseFloat(input_mass_element!.value)
    if (Number.isNaN(input_mass)) {
      alert("Input is not a Number")
      return;
    }

    table_body!.innerHTML += createTableRowHTML(input_mass)
}

function deleteRows() {
    const table = document.getElementById("ms1-datatable") as HTMLTableElement;
    const checkboxes = table!.querySelectorAll('tbody input[type="checkbox"]');

    checkboxes.forEach((checkbox) => {
        if ((checkbox as HTMLInputElement).checked) {
          const row = (checkbox.parentNode!.parentNode as HTMLTableRowElement);
          table!.deleteRow(row!.rowIndex);
        }
      });

}

function selectRows() {
    const table = document.getElementById('ms1-datatable') as HTMLTableElement;
    const selectAllCheckbox = table.querySelector('thead input[type="checkbox"]');
    const checkboxes = table.querySelectorAll('tbody input[type="checkbox"]');
  
    selectAllCheckbox!.addEventListener('change', (event) => {
      checkboxes.forEach((checkbox) => {
        (checkbox as HTMLInputElement).checked = (event.target as HTMLInputElement).checked;
      });
    });
}

export function createBottomRow(index: number, _names_for_img: string[], _smiles_for_img: string[]) {
  const row = document.createElement("tr");
  row.className = "hidden-row";
  row.id = "hidden-row-" + index.toString();

  return row;
}

function toggleRows(): void {
  const table: HTMLTableElement | null = document.getElementById('ms1-datatable') as HTMLTableElement;
  if (!table) {
      console.error('Table not found.');
      return;
  }

  const rows: NodeListOf<HTMLTableRowElement> = table.querySelectorAll('tbody tr.data');

  rows.forEach((row: HTMLTableRowElement) => {
      // Assuming the 3rd column is the one to check (index 2 since index starts at 0)
      const cell: HTMLTableCellElement = row.children[3] as HTMLTableCellElement;

      if (cell.textContent) {
        
          if (cell.textContent === '') {
              row.style.display = 'none';
          } else {
              row.style.display = 'table-row';
          }
      }
      else {
        if (row.style.display === "none") {
          row.style.display = '';
        } else {
          row.style.display = 'none';
        } 
        
      }
  });
}

