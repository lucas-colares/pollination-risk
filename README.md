# Food Plants in Brazil: Origin, Economic Value of Pollination and Pollinator Shortage Risk

This repository contains the data and R code used in the paper:

> **Food plants in Brazil: origin, economic value of pollination and pollinator shortage risk**  
> Oliveira et al. (2024)

The study provides an in‐depth analysis of 199 food crops in Brazil, focusing on their origins (native vs. exotic), pollinator dependence, and the economic value of pollination. It also examines the risk of pollinator shortages across different Brazilian regions and biomes. Key biological highlights include:
  
- **Pollinator interactions:** The Western honey bee (*Apis mellifera*) interacts extensively with both native and exotic crops, while vertebrate pollinators are found exclusively on native crops.
- **Pollinator dependence:** About 71.4% of native food crops are highly dependent on pollinators, in contrast to 30.2% of exotic crops.
- **Cultivation patterns:** Approximately 81.5% of Brazil’s agricultural area is devoted to exotic food crops, with soybean accounting for 46% of this area.
- **Regional risk:** Native crops, especially in the Northeast and Southeast, face a higher risk of yield loss from pollinator shortages.
- **Biome vulnerability:** Among Brazilian biomes, the Atlantic Forest is highlighted as particularly vulnerable to pollinator decline.

This research is especially significant for understanding how the expansion of monocultures and the heavy reliance on exotic crops might impact food security in Brazil—a megadiverse country and one of the world’s largest agricultural producers.

---

## Repository Structure

Below is an overview of the repository’s folder structure:

```
.
├── datasets
│   ├── all_spp.csv                # Lists all crop species with their taxonomic and identifier information.
│   ├── category.csv          # Defines categorical groupings (e.g., native vs. exotic) for the species.
│   ├── crops_info.csv          # Provides basic metadata and descriptive details for each crop.
│   ├── final_dataset_crops.csv          # The merged dataset combining crop, species, and pollinator dependency data for analysis.
│   ├── pol_dependency.csv          # Contains the pollinator dependency levels for each crop.
│   └── polinizadores.csv      # Lists the pollinator species data used in the study.
│
├── scripts
│   └── PAM - 2021           # Contains data from the Brazilian Agricultural Production (PAM) survey for 2021, detailing crop production statistics at the municipal level.
│       ├── area_colhida.csv    # Harvested area (in hectares) for different crops across Brazilian municipalities.
│       ├── area_plantada.csv    # Planted area (in hectares) for various crops in 2021.
│       ├── quant.csv    # Crop production quantity (in tons) reported at the municipal level.
│       ├── rend.csv    # Crop yield (kg per hectare) calculated for each crop.
│       └── valor.csv      # Economic value (in Brazilian Reais) of crop production.
│
├── PolEco.Rproj             # RStudio project file for organizing and managing the analysis environment, ensuring reproducibility and streamlined access to scripts and data.
│
└── README.md              # This file
```

---

## How to Reproduce the Analysis

1. **Data Preparation:**  
   The raw datasets in the `data/raw/` folder need to be preprocessed using the scripts found in `R/scripts/`. Running these scripts will generate cleaned datasets that are saved in `data/processed/`.

2. **Analysis Execution:**  
   Open the R Markdown files in the `R/analysis/` folder to follow the step-by-step analysis. These documents explain the statistical methods, models, and visualizations used in the paper.

3. **Figures:**  
   All visual outputs (e.g., maps, charts) will be automatically generated and stored in the `figures/` folder when you run the analysis scripts.

---

## Citation

If you use or modify any part of this repository in your work, please cite the original paper:

Oliveira, W., Colares, L.F., Porto, R.G., Viana, B.F., Tabarelli, M., & Lopes, A.V. (2024). *Food plants in Brazil: origin, economic value of pollination and pollinator shortage risk*. Science of The Total Environment, 169147. doi:10.1016/j.scitotenv.2023.169147

---

## License

This repository is released under the MIT license.
