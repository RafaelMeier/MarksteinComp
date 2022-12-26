# MarksteinComp

The script computes the Markstein number based on the strain rate of counter-flow flames.

The computation method follows https://www.sciencedirect.com/science/article/pii/S0010218002003681

If this code were used for any academic purpose, please reference it as

R. Meier, T. C. Cardoso, A. M. Oliveira, "Numerical Study of Triple Flames in Partially Premixed Methane and Hydrogen-Air Mixtures", _In: Proceedings of 19th Brazilian Congress of Thermal Sciences and Engineering_, Bento Gonçalves, Brazil, 2022.


## Libraries list

- Cantera >=2.5
- Pandas

## Pipeline

1. premixedCounterflow.py
2. postProcessingFlameStretch.py
3. postProcessingMarksteinNumber.py

#### OBS:
Some flames data for H2 and CH4 are provided in FlamesProps.py


