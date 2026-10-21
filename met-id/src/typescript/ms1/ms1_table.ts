
window.addEventListener("DOMContentLoaded", () => {
    document.getElementById("append-to-ms1-table")?.addEventListener("click", () => appendToMSTable());
    document.getElementById("ms1-table-delete-row")?.addEventListener("click", () => deleteRows())
    document.getElementById("ms1-table-head-check")?.addEventListener("click", () => selectRows())

    document.getElementById("toggle-show-identified")?.addEventListener("click", () => toggleRows())
    document.getElementById("toggle-show-significant")?.addEventListener("click", () => toggleRowsSignificant())

    const ionMobilitySwitch = document.getElementById(
        "toggle-ion-mobility"
    ) as HTMLInputElement | null;

    ionMobilitySwitch?.addEventListener(
        "change",
        updateIonMobilityVisibility
    );

    updateIonMobilityVisibility();
});

function updateIonMobilityVisibility(): void {
    const toggle = document.getElementById(
        "toggle-ion-mobility"
    ) as HTMLInputElement | null;

    if (!toggle) {
        return;
    }

    document.body.classList.toggle(
        "ion-mobility-off",
        !toggle.checked
    );
}



let table_body = document.getElementById("ms1-table-body");

export function createTableRowHTML(
    input_mz: number,
    input_ccs: number,
    include_ccs: boolean
): string {
    const ccsValue = include_ccs ? String(input_ccs) : "";

    return `
        <tr>
            <td><input type="checkbox" /></td>
            <td class="clickable-cell">${input_mz}</td>
            <td class="clickable-cell"></td>

            <td class="clickable-cell ion-mobility-column">
                ${ccsValue}
            </td>

            <td class="clickable-cell"></td>
            <td class="clickable-cell"></td>
            <td class="clickable-cell"></td>

            <td class="clickable-cell ion-mobility-column"></td>

            <td class="clickable-cell"></td>
            <td class="clickable-cell"></td>

            <td class="clickable-cell ion-mobility-column"></td>

            <td class="clickable-cell">
                <span class="unavailable"></span>
            </td>

            <td class="clickable-cell"></td>
            <td class="clickable-cell ion-mobility-column"></td>
            <td class="clickable-cell ion-mobility-column"></td>
        </tr>
    `;
}

function appendToMSTable() {
    let input_mass_element = document.getElementById("add-to-ms1-table-input") as HTMLInputElement;
    let input_ccs_element = document.getElementById("add-ccs-to-ms1-table-input") as HTMLInputElement;
    let input_mass = parseFloat(input_mass_element!.value)
    let input_ccs = parseFloat(input_ccs_element!.value)
    if (Number.isNaN(input_mass)) {
      alert("m/z is not a Number")
      return;
    }
    if (Number.isNaN(input_ccs)) {
      input_ccs = 0;
    }

    table_body!.innerHTML += createTableRowHTML(input_mass, input_ccs, true)
    //updateIonMobilityVisibility();
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
      // Assuming the 4th column is the one to check (index 2 since index starts at 0)
      const cell: HTMLTableCellElement = row.children[4] as HTMLTableCellElement;

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

function toggleRowsSignificant(): void {
    const table: HTMLTableElement | null = document.getElementById('ms1-datatable') as HTMLTableElement;
    if (!table) {
        console.error('Table not found.');
        return;
    }

    const rows: NodeListOf<HTMLTableRowElement> = table.querySelectorAll('tbody tr.data');

    rows.forEach((row: HTMLTableRowElement) => {
      // Assuming the 4th column is the one to check (index 2 since index starts at 0)
      const cell: HTMLTableCellElement = row.children[13] as HTMLTableCellElement;

      if (parseFloat(cell.textContent!) < 0.05) {
        
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

