# Feather Generator


A generative approach to creating down feathers



## Overall Process 

1. Draw a Curve in Rhino


2. Run the exportPath.py in python rhinoscript editor and select the curve


3. Open feather_v2.py in your IDE and go to the Hyperparameters area in the Main() function


4. Edit the hyper parameters (see documentation "Documentation.pdf") as you see fit and run the script


5. Run the importPaths.py in python rhinoscript editor. This will generate your output


6. Select output meshes and export as one .obj file. Wrap remesh file in z-brush or magicx for printing



## Built in relationships

* The further a barb is along its stem and the less barbs that barb will have (lower density)

* The further a barb is along its stem the thinner its initial radius

* The smaller the inital radius the shorter the barb is likely to be 