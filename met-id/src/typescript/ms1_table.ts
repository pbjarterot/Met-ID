
window.addEventListener("DOMContentLoaded", () => {
    document.getElementById("append-to-ms1-table")?.addEventListener("click", () => appendToMSTable());
    document.getElementById("ms1-table-delete-row")?.addEventListener("click", () => deleteRows())
    document.getElementById("ms1-table-head-check")?.addEventListener("click", () => selectRows())


});

function appendToMSTable() {
    console.log("button clicked");
    let input_mass_element = document.getElementById("add-to-ms1-table-input") as HTMLInputElement;
    let input_mass = input_mass_element!.value

    let table_body = document.getElementById("ms1-table-body");

    table_body!.innerHTML += `
                        <tr>
                          <td>
                            <input type="checkbox" />
                          </td>
                          <td>${input_mass}</td>
                          <td></td>
                          <td></td>
                          <td></td>
                          <td></td>
                          <td></td>
                          <td></td>
                        </tr>
                        `
}

function deleteRows() {
    console.log("trying to delete")
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
    /*
    checkboxes.forEach((checkbox) => {
      checkbox.addEventListener('change', () => {
        selectAllCheckbox!.checked = Array.from(checkboxes).every((checkbox) => (checkbox as HTMLInputElement).checked);
      });
    });
    */
  }
