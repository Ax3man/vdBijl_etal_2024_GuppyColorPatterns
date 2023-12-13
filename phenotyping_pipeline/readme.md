---
editor_options: 
  markdown: 
    wrap: 72
---

# The phenotyping pipeline

This folder contains the scripts used to ingest, annotate and transform the guppy photos. The code used to actually train the deep learning models used in the pipeline can be found within the `phenotyping_models/` folder. The rest of the files here contain the code to build and run the pipeline itself.

`new_photo_intake.R` was ran to detect newly added photos, run all the pipeline components on those new photos, add them to the database, and generate data files.

`generate_consensus_shape.R`, `generate_warped_images.R` and `process_full_landmarks.R` were run to programmatically generate the consensus shape, the images warped to the consensus shape, and the sets of "full landmarks" (along the whole outline of the fish), respectively.

These files are used by the scripts above and contain functions and tools: 

- `tools.R`
- `index_photos.R`
- `fish_extraction.R`
- `auto_place_landmarks.R`
- `carotenoid_extraction.R`
- `melanin_extraction_v2.R`
- `place_manual_landmarks.R`
