# 💸 How to price voltage stability ancillary services for adequate profits

This repository contains Julia files for pricing **Voltage Stability Ancillary Services via Primal-Dual Formulation** in IBG-dominated power systems, using the **JUMP** package. Three pricing approaches have been examined on addressing non-convexities from Unit Commitment (UC):

1. **Restricted pricing method** – Based on duality theory  
2. **Dispatchable pricing method** – Based on duality theory
3. **Primal-dual formulation** – Based on duality theory
   
---

# 📖 Relevant works

The main models and methodologies are in the listed papers here.
1. Voltage stability constraints refer to:
- [Zhongda Chu, and Fei Teng. "Voltage Stability Constrained Unit Commitment in Power Systems With High Penetration of Inverter-Based Generators." IEEE Transactions on Power Systems.](https://ieeexplore.ieee.org/abstract/document/9786660)
2. Pricing approaches refer to:
- [Peng Wang, and Luis Badesa. "Shadow Pricing of Static Voltage Stability Services within Unit Commitment for Inverter-Dominated Power Systems." arxiv](https://arxiv.org/abs/2607.05209)
- [Peng Wang, and Luis Badesa. "Pricing Short-Circuit Current via a Primal-Dual Formulation for Preserving Integrality Constraints." arxiv](https://arxiv.org/abs/2510.05293)
- [Gribik PR, Hogan WW, and Pope SL. "Market-Clearing Electricity Prices and Energy Uplift." Cambridge, MA](https://www.academia.edu/download/45461778/hogan2.pdf)
3. Data used in this work and relevant work/data refer to:
- [Wang, Peng, and Luis Badesa. "Imperfect Competition in Markets for Short-Circuit Current Services." Sustainable Energy, Grids and Networks (2026): 102423.](https://arxiv.org/pdf/2508.09425)

---

# 🔧 Guidance on how to use the code

The work is mainly made of two parts:

1. **Modelling voltage stability constraints.**
2. **Modelling pricing schemes.**

- The main function including modelling voltage stability constraints and pricing approaches is `Main.jl`.
- For the code of voltage stability constraints modelling, please refer to the files named `admittance_matrix_calculation.jl`, `dataset_gene.jl` and `offline_trainning.jl`.
  1. `admittance_matrix_calculation.jl` calculates the impedance of transmission lines of the system.
  2. `dataset_gene.jl` generates the data for classification, i.e., the offline trainning process. The subfunction `admittance_matrix_calculation.jl` is called here to obtain the transmission line admittance matrix which is combined with the generators' admittance matrix. In `dataset_gene.jl`, the actual representation of nonlinear termes mentioned in the paper are included. The remainder of code is generating all possible operating states of generators.
  3. `offline_trainning.jl` is the trainning process, with inputting parameters from above subfunctions.
- For the code of pricing market services, please refer to the files named `pricing_VS_restricted.jl`, `pricing_VS_dispatchable.jl` and `pricing_VS_PD_fixprimal_multi_period.jl`.
  4. `reve_pro_calculation_dispatchable.jl` and `reve_pro_calculation_restricted_approx.jl` are the subfunctions for computing revenues and profits of generating units.

----
# ✒️ Citation
If you find something helpful or use this code for your own work, please cite this paper:
<ol>
      Wang, Peng and Luis Badesa. "A Primal-Dual Formulation for Pricing Static Voltage Stability Services within a Unit Commitment Model." arXiv preprint arXiv:- (2026).
</ol>
      <br>
      
<ol> 

```bibtex
@misc{-,
  title        = {A Primal-Dual Formulation for Pricing Static Voltage Stability Services within a Unit Commitment Model},
  author       = {Peng, Wang and Luis, Badesa},
  year         = {2026},
  eprint       = {-},
  archivePrefix= {-},
  primaryClass = {-},
  url          = {-}
}
```

---

# 📃 Funding

This work was supported by MICIU/AEI/10.13039/501100011033 and ERDF/EU under grant PID2023-150401OA-C22, as well as by the Madrid Government (Comunidad de Madrid-Spain) under the Multiannual Agreement 2023-2026 with Universidad Politécnica de Madrid, `Line A - Emerging PIs' (grant number: 24-DWGG5L-33-SMHGZ1). The work of Peng Wang was also supported by China Scholarship Council under grant 202408500065.

<figure style="display:inline-block; margin:10px; text-align:center;">
   <img src="./logos/MICIU+Cofinanciado+AEI.jpg" style="width:600px; height:140px; object-fit:contain; display:block;">
</figure>
<figure style="display:inline-block; margin:10px; text-align:center;">
   <img src="./logos/Logo_CM.png" style="width:200px; height:160px; object-fit:contain; display:block;">
</figure>
