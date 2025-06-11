import * as d3 from 'd3';
import { invoke } from '@tauri-apps/api/core';
import { convertTableToCSV, importFromCSV } from './ms1_io';

function mass_error_slidein() {
    const slideInContent = document.getElementById("slide-in-bottom-content");
    if (slideInContent) {
      slideInContent.innerHTML += `
        <div class="manual-mass-error-top-div">
          <div class="mass-error-manual-input-div-theoretical">
            <p class="mass-error-manual-text">Theoretical m/z</p>
            <input class="mass-error-manual-input" id="mass-error-manual-input-theoretical" placeholder="0.0"></input>
          </div>
          <div class="mass-error-manual-input-div-observed">
            <p class="mass-error-manual-text">Observed m/z</p>
            <input class="mass-error-manual-input" id="mass-error-manual-input-observed" placeholder="0.0"></input>
          </div>
          <!-- 
          <div class="mass-error-manual-input-div-lock">
            <p class="mass-error-manual-text">Lockmass m/z</p>
            <input class="mass-error-manual-input" id="mass-error-manual-input-lock" placeholder="0.0"></input>
          </div>
          -->
          <div class="ms1-manual-error-button-container">
            <button class="ms1-manual-error-button" id="ms1-manual-error-add-button">
              <span class="ms1-manual-error-button-text">Add to list</span> 
              <span class="ms1-manual-error-button-icon">
              <ion-icon name="add-circle-outline"></ion-icon>
              </span>
            </button>
          </div>
          <!-- 
          <div class="ms1-manual-error-button-container">
            <button class="ms1-manual-error-button" id="ms1-manual-error-add-button">
              <span class="ms1-manual-error-button-text">Use List</span> 
              <span class="ms1-manual-error-button-icon">
              <ion-icon name="arrow-forward-circle-outline"></ion-icon>
              </span>
            </button>
          </div>
          -->
          <div class="ms1-manual-error-button-container">
            <button class="ms1-manual-error-button" id="ms1-manual-error-import-button">
              <span class="ms1-manual-error-button-text">Import List</span> 
              <span class="ms1-manual-error-button-icon">
              <ion-icon name="arrow-up-circle-outline"></ion-icon>
              </span>
            </button>
          </div>
          <div class="ms1-manual-error-button-container">
            <button class="ms1-manual-error-button" id="ms1-manual-error-export-button">
              <span class="ms1-manual-error-button-text">Export List</span> 
              <span class="ms1-manual-error-button-icon">
              <ion-icon name="arrow-down-circle-outline"></ion-icon>
              </span>
            </button>
          </div>
        </div>
      </div>
        <div class="manual-mass-error-table-div">
            <table class="ms1-mass-error-table" id="ms1-mass-error-table">
                <thead>
                    <tr>
                        <th>m/z</th>
                        <th>PPM diff</th>
                        <th></th>
                    </tr>
                </thead>
                <tbody id="ms1-mass-error-table-body">
                    <!-- Add rows as needed -->
                </tbody>
            </table>
        </div>
        <div class="manual-mass-error-plot-div" id="manual-mass-error-plot-div">
        </div>
      `;
    }


    const table = document.querySelector('.ms1-mass-error-table');
    table?.addEventListener('click', function(event: Event) {
        // Check if the clicked element is the trash can icon
        if ((event.target as HTMLElement)?.classList.contains('remove-mass-error-row')) {
        // Get the row element associated with the clicked trash can icon
        const row = (event.target as HTMLElement).closest('tr');
    
        if (row && row.parentNode) {
            // Remove the row from the table
            row.parentNode.removeChild(row);
        }
        }
    });
    const plotData = new PlotData("#manual-mass-error-plot-div");

    const addErrorPointButton = document.getElementById("ms1-manual-error-add-button");
    const importErrorPointsButton = document.getElementById("ms1-manual-error-import-button");
    const exportErrorPointsButton = document.getElementById("ms1-manual-error-export-button");
    const theoretical_err = document.getElementById("mass-error-manual-input-theoretical") as HTMLInputElement;
    const observed_err = document.getElementById("mass-error-manual-input-observed") as HTMLInputElement;
    addErrorPointButton?.addEventListener("click", () => {
      plotData.addPoint(parseFloat(theoretical_err.value), parseFloat(observed_err.value), true);
    });

    importErrorPointsButton?.addEventListener("click", async () => {
      let data: string[][] = await importFromCSV();
      add_to_error_table(data, plotData);
    })

    exportErrorPointsButton?.addEventListener("click", () => {
      convertTableToCSV("ms1-mass-error-table");
    })
  }
function parse_error_number(num: string) {
  const parsedNumber: number = parseFloat(num);
  if (!isNaN(parsedNumber)) {
    return parsedNumber;
} else {
    return 0;
}
}
  

function add_to_error_table(file_contents: string[][], plotData: PlotData) {
    const tbody = document.getElementById("ms1-mass-error-table-body")
        tbody!.innerHTML = ""

        for (let i = 0; i < file_contents.length; i++) {
            let dat: Datum = {x: parse_error_number(file_contents[i][0]), y: parse_error_number(file_contents[i][1])}
            plotData.addPoint(dat.x, dat.y, false);
        }
}
  


const slideInDiv = document.getElementById("slide-in-bottom") as HTMLElement;
const resizeHandle = document.querySelector(".resize-handle") as HTMLElement;
const slideinButton = document.getElementById("ms1-error-button-manual");
const plotDiv = document.querySelector('.manual-mass-error-plot-div') as HTMLElement;

slideinButton?.addEventListener("click", () => {
  slideInDiv.classList.toggle("slide-in-bottom--active");
  if (slideInDiv.classList.contains("slide-in-bottom--active")) {
    mass_error_slidein();
  } else {
    const slideInContent = document.getElementById("slide-in-bottom-content");
    if (slideInContent) {
      slideInContent.innerHTML = "";
    }
  }
});
  
resizeHandle.addEventListener("dragstart", (event: DragEvent) => {
  event.dataTransfer!.setDragImage(new Image(), 0, 0); // Disables default drag image
  slideInDiv.classList.add("resizing");
});

resizeHandle.addEventListener("drag", (event: DragEvent) => {
  const mouseY = event.clientY;
  const newHeight = window.innerHeight - mouseY;
  slideInDiv.style.height = newHeight + "px";
});

resizeHandle.addEventListener("dragend", () => {
  slideInDiv.classList.remove("resizing");
});

plotDiv?.addEventListener("click", () => {
  slideInDiv.classList.toggle("plot-div-expanded");
});

slideInDiv.addEventListener("transitionend", () => {
  if (slideInDiv.classList.contains("slide-in-bottom--active")) {
    slideInDiv.style.transition = "none"; // Disable transition after sliding in
  }
});






interface Datum {
  x: number;
  y: number;
}

class PlotData {
  private data: Datum[] = [];
  private parentDiv: HTMLDivElement;
  private parentParentDiv: HTMLDivElement;
  private svg: d3.Selection<SVGSVGElement, unknown, HTMLElement, any>;
  private xScale!: d3.ScaleLinear<number, number>; // Use "!" to assert definite assignment
  private yScale!: d3.ScaleLinear<number, number>; // Use "!" to assert definite assignment
  private margin: {left: number, right: number, top: number, bottom: number};
  private xMin: number = 0;
  private xMax: number = 100;
  private yMin: number = 0;
  private yMax: number = 100;
  private plotWidth: number;
  private plotHeight: number;


  constructor(parentDivSelector: string) {
    this.parentDiv = document.querySelector(parentDivSelector) as HTMLDivElement;
    this.parentParentDiv = this.parentDiv.parentElement as HTMLDivElement;

    this.margin = { top: 20, right: 20, bottom: 30, left: 40 };
  
    // Create the SVG element
    const svgElement = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    svgElement.setAttribute("class", "scatterplot");
  
    // Append the SVG element to the parent div
    this.parentDiv.appendChild(svgElement);
  
    // Select the SVG element using d3.select
    // Adjust the parent type in the selection to match with HTMLElement
    this.svg = d3.select<SVGSVGElement, unknown>(svgElement as SVGSVGElement) as unknown as d3.Selection<SVGSVGElement, unknown, HTMLElement, any>;
  
    // Modify the type of the data argument in the data() method
    this.svg.data<unknown[]>([this.data]);

    this.plotWidth = document.getElementById("manual-mass-error-plot-div")!.clientWidth; - this.margin.left - this.margin.right;
    this.plotHeight = document.getElementById("manual-mass-error-plot-div")!.clientHeight; - this.margin.top - this.margin.bottom;
  
    // Call updateSize initially to set up the plot
    this.updateSize(200);
    this.updateSize = this.updateSize.bind(this);
    this.updateScatterplot = this.updateScatterplot.bind(this);
    this.addPoint = this.addPoint.bind(this);
    this.ppmDifference = this.ppmDifference.bind(this);
    this.updateScale = this.updateScale.bind(this);
  
    // Attach event listeners to the parent's parent div
    this.parentParentDiv.addEventListener('mouseenter', this.onMouseEnter.bind(this));
    this.parentParentDiv.addEventListener('mouseleave', this.onMouseLeave.bind(this));
  }

  private onMouseEnter() {
    this.updateSize(400);
  }
  private onMouseLeave() {
    this.updateSize(200);
  }

  public updateSize(height: number): void {
    const width = document.getElementById("manual-mass-error-plot-div")!.clientWidth;
  
    // Set the width and height of the SVG element
    this.svg
      .transition()
      .duration(200)
      .attr("width", width)
      .attr("height", height);
  
    // Define margins
    this.margin = { top: 20, right: 20, bottom: 40, left: 40 };
    this.plotWidth = width - this.margin.left - this.margin.right;
    this.plotHeight = height - this.margin.top - this.margin.bottom;
  
    // Set the width and height of the plot area within the SVG
    let plotArea = this.svg.select<SVGGElement>('g.plot-area');
    if (plotArea.empty()) {
      plotArea = this.svg.append('g').attr('class', 'plot-area');
    }
    plotArea.attr("transform", `translate(${this.margin.left}, ${this.margin.top})`);
  
    // Update scales based on data
    if (this.data.length > 0) {
      this.xMin = d3.min(this.data, (d) => d.x) || 0;
      this.xMax = d3.max(this.data, (d) => d.x) || 0;
      this.yMin = d3.min(this.data, (d) => d.y) || 0;
      this.yMax = d3.max(this.data, (d) => d.y) || 0;
    }
  
    this.xScale = d3.scaleLinear().domain([this.xMin - 10, this.xMax + 10]).range([0, this.plotWidth]);
    this.yScale = d3.scaleLinear().domain([this.yMin - 1, this.yMax + 1]).range([this.plotHeight, 0]);
  
    // Remove existing axes
    this.svg.select<SVGSVGElement>("g.x-axis").remove();
    this.svg.select<SVGSVGElement>("g.y-axis").remove();
  
    // Create x-axis
    const xAxis = d3.axisBottom(this.xScale)
      .tickSizeOuter(0)
      .tickSizeInner(-this.plotHeight)
      .tickPadding(10)
      .tickFormat(d3.format("d"));
  
    const xAxisGroup = this.svg
      .append("g")
      .attr("class", "x-axis")
      .attr("transform", `translate(${this.margin.left}, ${this.margin.top + this.plotHeight})`)
      .call(xAxis);
  
    xAxisGroup
      .selectAll(".tick text")
      .style("fill", "white")
      .style("font-size", "12px");
  
    // Create y-axis
    const yAxis = d3.axisLeft(this.yScale)
      .tickSizeOuter(0)
      .tickSizeInner(-this.plotWidth)
      .tickPadding(10)
      .tickFormat(d3.format("d"));
  
    const yAxisGroup = this.svg
      .append("g")
      .attr("class", "y-axis")
      .attr("transform", `translate(${this.margin.left}, ${this.margin.top})`)
      .call(yAxis);
  
    yAxisGroup
      .selectAll(".tick text")
      .style("fill", "white")
      .style("font-size", "12px");
  
    // Call updateScatterplot to update the plot
    this.updateScatterplot();
    this.updateScale();
    this.fitEquation();
  }

  public updateScatterplot(): void {
    // Bind data
    const circles = this.svg.selectAll<SVGCircleElement, Datum>("circle").data(this.data);
  
    // Enter new data
    circles.enter()
      .append("circle")
      .attr("r", 3)
      .attr("fill", "white")
      .merge(circles) // Merge enter and existing circles
      .attr("cx", (d) => this.xScale(d.x) + this.margin.left)
      .attr("cy", (d) => this.yScale(d.y) + this.margin.top);
  
    // Remove old data
    circles.exit().remove();
  
    // Update the x-axis
    const xAxis = d3.axisBottom(this.xScale);
    this.svg.select<SVGSVGElement>("g.x-axis").call(xAxis)
      .selectAll(".tick text")
      .style("fill", "white")
      .style("font-size", "12px");
  
    // Update the y-axis
    const yAxis = d3.axisLeft(this.yScale);
    this.svg.select<SVGSVGElement>("g.y-axis").call(yAxis)
      .selectAll(".tick text")
      .style("fill", "white")
      .style("font-size", "12px");

    //this.fitEquation();
  }
  
  private updateScale(): void {
  
    if (this.data.length > 0) {
      this.xMin = d3.min(this.data, (d) => d.x as number) as number || 0;
      this.xMax = d3.max(this.data, (d) => d.x as number) as number || 0;
      this.yMin = d3.min(this.data, (d) => d.y as number) as number || 0;
      this.yMax = d3.max(this.data, (d) => d.y as number) as number || 0;
    }
  
    this.xScale.domain([this.xMin as number -10, this.xMax as number +10]).range([0, this.plotWidth]);;
    this.yScale.domain([this.yMin as number -1, this.yMax as number +1]).range([this.plotHeight, 0]);;
  
    // Update the x-axis
    const xAxis = d3.axisBottom(this.xScale);
    this.svg.select<SVGSVGElement>("g.x-axis").call(xAxis);
  
    // Update the y-axis
    const yAxis = d3.axisLeft(this.yScale);
    this.svg.select<SVGSVGElement>("g.y-axis").call(yAxis);
  }

  public addPoint(xValue: number, yValue: number, calcppm: boolean): void {
    let newDataPoint: Datum;
    if (calcppm) {
      let newYvalue = this.ppmDifference(xValue, yValue);
      newDataPoint = { x: xValue, y: newYvalue };
    } else {
      newDataPoint = { x: xValue, y: yValue };
    }
    this.data.push(newDataPoint);
  
    this.updateScale();
    this.updateScatterplot();
    this.addToTable(newDataPoint);
    this.fitEquation();
  }
  
  public addToTable(dataPoint: Datum): void {
    const table = document.getElementById("ms1-mass-error-table") as HTMLTableElement;
    if (table) {
      const row = table.insertRow();
      const xCell = row.insertCell();
      const yCell = row.insertCell();
      const deleteCell = row.insertCell();
  
      xCell.innerHTML = dataPoint.x.toString();
      yCell.innerHTML = dataPoint.y.toString();
  
      const deleteButton = document.createElement("button");
      const trashIcon = document.createElement("ion-icon") as HTMLElement;
      trashIcon.setAttribute("name", "trash-outline"); // Set the name of the Ionicon (trash outline)
      trashIcon.setAttribute("slot", "icon-only"); // Use the icon-only slot for better styling
      deleteButton.appendChild(trashIcon);
      deleteButton.addEventListener("click", () => {
        this.deletePoint(dataPoint);
        table.deleteRow(row.rowIndex);
        this.updateScale()
        this.fitEquation();
      });
  
      deleteCell.appendChild(deleteButton);
    }
  }
  
  private deletePoint(dataPoint: Datum): void {
    const index = this.data.findIndex((d) => d === dataPoint);
    if (index !== -1) {
      this.data.splice(index, 1);
      this.updateScale();
      this.updateScatterplot();
      this.fitEquation();
    }
  }

  private ppmDifference(theoretical: number, observed: number): number {
    if (theoretical === 0 || observed === 0) {
      throw new Error("value1 cannot be 0");
    }
    const difference =  observed - theoretical;
    const ppmDifference = (difference / theoretical) * 1_000_000;
    return ppmDifference;
  }

  public async fitEquation(): Promise<void> {
    let mass_error_input = document.getElementById("ms1-error-input-text") as HTMLInputElement;


      if (this.data.length > 0) {
        const tuples: [number, number][] = this.data.map(datum => [datum.x, datum.y])
      let linear_coeffs: number[] = await invoke("mass_error_regression", {data: tuples})
      
      mass_error_input!.value = linear_coeffs.toString();

      // Update the plot with the linear equation
      this.updatePlotWithQuadraticEquation(linear_coeffs);
      }
      
    
  }

  private updatePlotWithQuadraticEquation(equation: number[]): void {
    // Remove any existing line or curve
    this.svg.select('.equation-line').remove();
    this.svg.select('.equation-curve').remove();

    const a = equation[0];
    const b = equation[1];
    const c = equation[2];
    

    // Generate a curve using the quadratic equation
    const xValues = d3.range(this.xScale.domain()[0], this.xScale.domain()[1], 0.1);
    const curveData = xValues.map((x) => ({ x, y: a * x * x + b * x + c }));

    // Create a curve generator
    const curve = d3.line<{ x: number; y: number }>()
      .x((d) => this.xScale(d.x))
      .y((d) => this.yScale(d.y));

    // Append the curve to the plot area
    this.svg.select<SVGSVGElement>('g.plot-area')
      .append('path')
      .attr('class', 'equation-curve')
      .datum(curveData)
      .attr('d', curve)
      .style('stroke', 'red')
      .style('stroke-width', '2')
      .style('fill', 'none');
  }
}






  



  
  





