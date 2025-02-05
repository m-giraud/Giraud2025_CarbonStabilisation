# CPlantBox
## Coupled models
Coupling of the [CPlantBox model](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox) and its [dumux-rosi module](https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi) (based on [DuMux]([http://cplantbox.com](https://dumux.org/)))


## Setup the files

1) Download the files in a linux environment.
2) make sure you have all the requirements installed by running the file "checkRequirements.py"
## Recreate the main results
run the 9 simulations presented in the paper (baseline, earlry dry spell, late dry spell for biokinetic parameter set 5 44 and 61), by running:
```
cd dumux38TraiRhizo/dumux/dumux-rosi/experimental/fixedPointIter2/scripts 

python3 mainTraiRhizo.py 10 25 5 baseline  
python3 mainTraiRhizo.py 10 25 5 earlyDry
python3 mainTraiRhizo.py 10 25 5 lateDry
python3 mainTraiRhizo.py 10 25 44 baseline  
python3 mainTraiRhizo.py 10 25 44 earlyDry  
python3 mainTraiRhizo.py 10 25 44 lateDry 
python3 mainTraiRhizo.py 10 25 61 baseline
python3 mainTraiRhizo.py 10 25 61 earlyDry
python3 mainTraiRhizo.py 10 25 61 lateDry
```
