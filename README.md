# MarksteinComp

The MarksteinComp scripts compute the Markstein number based on the strain rate of premixed counter-flow flames. 

For weakly-stretched flames, the Markstein length is a phenomenological parameter that depends on properties such as the unburnt and burnt mixture temperature, local heat release, and the Lewis number. It characterizes a linear variation of the flame speed associated with the local flame surface topology regarding the local flame front curvature and/or strain rate.

The computation method follows https://www.sciencedirect.com/science/article/pii/S0010218002003681

If this code were used for any academic purpose reference it as

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


