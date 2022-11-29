//______________________________________Collagen measurenments around vessels_____________________//
// Markers: COL1A1:Collagen, CD31:Endothelial Cells, PDGFR-B: Pericytes, GFAP:tumor cells and astrocytes
//rename image and split channels

rename("sample");
run("Stack to Images");

selectWindow("FITC");
rename("Collagen");

selectWindow("Aqua");
rename("Endothelial");

selectWindow("Cy5");
rename("Pericytes");

selectWindow("Red610");
rename("Tumor");


//___________________________________________Measuring Average GFAP___________________________________________//

//Measure the average signal intensity of the tumor (GFAP) marker in each tile
selectWindow("Tumor");
run("Red");
run("Set Measurements...", "area mean modal min integrated display redirect=None decimal=3");
run("Measure");
Table.rename("Results", "Results-Tumor");

//___________________________________________Detecting the Vessels___________________________________________//

selectWindow("Endothelial");

// Create a duplicate window in which the vessel area is detected, this window will be named "Endothelial-1", keep the original window for expoting data or further assessments
run("Duplicate...", " ");
run("Enhance Contrast", "saturated=5.35");
run("Gaussian Blur...", "sigma=5");

// CD31 signal that is above a ceratin threshold (called a) will be included as possibly being part of a vessel, any signal below this value is excluded.
//This threshold is set based on the mean, maximum, and standard deviation values of each image. 
//This lower threshold is set to a=mean+b*c, where b depend on mean and standard deviation values, and c depends on maximum  and standard deviation values.

b=7;
getStatistics(area, mean, min, max, std);

if (mean > 43) {
	if (std<50) {
		b=19;
	}
		if (std>50 && std<60) {
		b=16;
	}	
	if (std>60) {
		b=12;
	}
}

if (mean > 50) {
	if (std<50) {
		b=27;
	}
		if (std>50 && std<60) {
		b=22;
	}
		if (std>60) {
		b=16;
	}
}

if (mean > 80) {
	if (std<60) {
		b=32;
		}
	if (std>60) {
		b=25;
		}
}

if (max<300) {
		c=1;
	}
if (max>300 && max<1000) {
		c=1.5;
	}
if (max>1000) {
		c=2;
	}

if (mean < 43){
	c=1;
	if (std < 38) {
		b=9;}
	if (std < 30) {
		b=10;}
	if (std < 18) {
		b=12;}
	if (std < 12) {
		b=16;}
}

a=mean+c*b;

// After the lower threshold is found, a binary mask will be created.
setThreshold(a, 65535, "raw");
setOption("BlackBackground", true);
run("Convert to Mask");

// Only areas that are larger than 150 pixels will be included as individual vessels.
run("Analyze Particles...", "size=150-Infinity clear include add");
selectWindow("Endothelial");
roiManager("Show All");

//______________________________________Creating Rings Around the Vessels___________________________________//
// Four rings around each vessel are created
rings   = 4;
//The distance between these rings is set to 3.5
diameter= 3.5; 

//Rename the oroginal vessel selection to organize them
n = roiManager('count');
for (i = n-1; i > -1; i--) {
    roiManager('select', i);
    roiManager("Rename", "Start_" + (i + 1));
    }

//Set measurements to include shape and Ferets measurements of each vessel
run("Set Measurements...", "area shape feret's redirect=None decimal=3");

//Dependng on the original size of the detected vessels, make all selections smaller to include the endothelial cell regions in collagen and pericyte measurements.  

count = roiManager("count");
array = newArray(count);

for (i=0; i<array.length; i++) {
   roiManager("select", i);
   roiManager("Measure");
   mferet = getResult("MinFeret", 0);
   area = getResult("Area", 0);
   if ((mferet > 20) && area>800) {
   run("Enlarge...", "enlarge=-5");
   Roi.setName("Circle_" + (i + 1) + "_" + "1")
   roiManager("Add");
   } else if (mferet > 14 && area>200){
   run("Enlarge...", "enlarge=-3.5");
   Roi.setName("Circle_" + (i + 1) + "_" + "1");
   roiManager("Add");
   	} else {
   	run("Enlarge...", "enlarge=-1");
   Roi.setName("Circle_" + (i + 1) + "_" + "1");
   roiManager("Add");
   }
   run("Clear Results");
}

//Reset measurements
run("Set Measurements...", "area mean min display redirect=None decimal=3");

//Delete original vessel selections from roiManager
n = roiManager('count');
for (i = n-1; i > -1; i--) {
    roiManager('select', i);
    name = Roi.getName();
    a = substring(name,0,5);
    if (matches(a, "Start")) {
    roiManager("delete");
    }
}

//Add circles around new (smaller) ROIs 
n = roiManager('count');
for (j = 0; j < n; j++) {
	   roiManager("Select",  j);
   for (i = 0; i <rings; i++){
	   run("Enlarge...", "enlarge=" + diameter);
	   Roi.setName("Circle_" + (j + 1) + "_" + (i + 2));
	   roiManager("Add");
   }
}

//Arrange ROIs by name
roiManager("sort");

//Make rings into "donut" (regions between the rings) ROIs to perform the collagen signal and pericyte signal measurements in the surfaces between the rings
n = roiManager('count');
for (i = 0; i < n; i += 5) {
	tmp = newArray(i,i+1,i+2,i+3,i+4);
   for (j = tmp[0]; j <tmp[4]; j++){
	   roiManager("Select", newArray(j,j+1));
	   roiManager("XOR");
       Roi.setName((i/5) + "_Ring_" + (j+1));
	   roiManager("Add");
   }
}

//Delete original circles from ROI manager
n = roiManager('count');
for (i = n-1; i > -1; i--) {
    roiManager('select', i);
    name = Roi.getName();
    a = substring(name,0,6);
    if (matches(a, "Circle")) {
    roiManager("delete");
    }
}

//Measure Collagen signal
run("Set Measurements...", "area mean modal min integrated display redirect=None decimal=3");
roiManager("Show None");
roiManager("Show All");
selectWindow("Collagen");
roiManager("Measure");
Table.rename("Results", "Results-Collagen-in-Rings");

//______________________________________________Cell count _____________________________________________________//

//Counting the nuclei:

//ROI was cleared from previous selections.
 if (roiManager("count")==0){
	} else {
			roiManager("Show All");
			roiManager("Deselect");
			roiManager("Delete");
	}
// At times, CD31 channel bled through the DAPI channel. To correct for this issue, the CD31 signal was subtracted from the DAPI signal before cell counting.
selectWindow("Endothelial");
roiManager("Show None");
run("Enhance Contrast...", "saturated=2");
run("Gaussian Blur...", "sigma=4");
imageCalculator("Subtract create", "DAPI","Endothelial");

//Cells were counted with Stardist using the DAPI signal.
selectWindow("Result of DAPI");
run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'Result of DAPI', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'25.0', 'percentileTop':'100.0', 'probThresh':'0.75', 'nmsThresh':'0.45', 'outputType':'ROI Manager', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
n = roiManager('count');
print("The cell count is: \n");
print(n);

