

export function addSearchbarListener() {
  const searchbarDiv = document.getElementById("ms1-searchbar-input") as HTMLInputElement;

  if (searchbarDiv) {
    searchbarDiv.addEventListener("input", (event) => {
      if (event.target instanceof HTMLInputElement) {
        MS1Searchresults(searchbarDiv!.value);
      }
    })
  }
}

function MS1Searchresults(searchvalue: string) {
  const rows = document.querySelectorAll('table tr');

    rows.forEach(row => {
        // Check if 'foo' is not contained in the row
        if (!row!.textContent!.includes(searchvalue)) {
            // Hide the row
            (row as HTMLTableRowElement).style.display = 'none';
        } else {
          (row as HTMLTableRowElement).style.display = 'table-row';
        }
    });
}
