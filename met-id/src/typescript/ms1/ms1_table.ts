window.addEventListener("DOMContentLoaded", () => {
    document.getElementById("append-to-ms1-table")?.addEventListener("click", () => appendToMSTable());
    document.getElementById("ms1-table-delete-row")?.addEventListener("click", () => deleteRows())
    document.getElementById("ms1-table-head-check")?.addEventListener("click", () => selectRows())


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
    let input_mass = input_mass_element!.value

    table_body!.innerHTML += createTableRowHTML(parseFloat(input_mass))
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
