Large fitted Bayesian model objects are hosted separately on Zenodo because several exceed GitHub's file-size limits.

To reproduce the analyses without refitting the JAGS models, download all .rds files form the associated Zenodo upload and place them directly in this results/ folder.

The analysis scripts automatically load the corresponding fitted object from this folder when available. If a required object is missing, the relevant model will be refitted and saved here. 