###SDM - Ruppells vulture

setwd ("~/Downloads/Ruppell's vulture")

## load the required packages
library(biomod2)
library(raster)
library(rasterVis)
library(gridExtra)
library(reshape2)

Ruppells <- read.csv("/Users/peggymutheu/Downloads/Ruppell's vulture/Ruppellsvulture_cleaned.csv1/Cleaned_ruppells.csv", header = TRUE)
#Loading predictor variables
lst <- list.files(path="~/Downloads/Ruppell's vulture/Variables",pattern='asc$',all.files = TRUE, full.names = T) 
preds<-stack(lst)

plot(preds[[1]])



#######################################################################################################
## curent climatic variables
stk_current <- 
  raster::stack(
    c(preds),
    RAT = FALSE
  )
plot(stk_current)

# Rename the layer names
new_names <- c("aspect", "elevation", "Human influence index", "NDVI", "Slope", "bio13", "bio14", "bio18", "bio19", "bio3", "bio7")  # Specify the new names as a character vector
names(stk_current) <- new_names
print(names(stk_current))
plot(stk_current)
occurrence_data <- Ruppells
head (occurrence_data)
occurrence_data$Presence <- 1
occurrence_data
data <- occurrence_data
table (data$species)
plot(data, add= T)
spp_to_model <- unique(data$species)
class(data)
head(data)
bioclim_data <- stk_current


biomod2_wrapper <- function(sp){
  cat("\n> species : ", sp)
  
  ## get occurrences points
  sp_dat <- data[data$species == sp, ]
  
  ## formatting the data
  myBiomodData <- 
    BIOMOD_FormatingData ( 
      resp.name = sp,
      resp.var = sp_dat$Presence, 
      expl.var = bioclim_data,
      resp.xy = sp_dat[, c("decimalLongitude", "decimalLatitude")],
      PA.strategy = "random", 
      PA.nb.rep = 2, 
      PA.nb.absences = 10000
    )
  ## print formatting summary
  myBiomodData
  if(!exists(sp)) dir.create(sp)
  pdf(paste(sp, "/", sp ,"_data_formated.pdf", sep="" ))
  try(plot(myBiomodData))
  dev.off()
  
  myBiomodOptions <- BIOMOD_ModelingOptions (GLM = list(type = 'quadratic', 
                                                        interaction.level = 1),
                                             GBM = list(n.trees = 5000),
                                             GAM = list(type = 's_smoother',
                                                        interaction.level = 0),
                                             RF = list (n.tree = 750))
  
  myBiomodOptions
  
  ## model species
  myBiomodModelOut <- BIOMOD_Modeling(
    bm.format = myBiomodData,
    modeling.id = 'AllModels',
    models = c("GBM", "GLM", "RF", "GAM"),
    bm.options = myBiomodOptions,
    CV.strategy = 'random',
    CV.nb.rep = 5,
    CV.perc = 0.7,
    CV.do.full.models = TRUE,
    metric.eval = c('TSS', 'ROC'),
    var.import = 3
  )
  
  myBiomodModelOut
  
  ###Get evaluation scores & variables importance
  get_evaluations(myBiomodModelOut)
  
  ## build ensemble models
  myBiomodEM <- BIOMOD_EnsembleModeling(
    bm.mod =  myBiomodModelOut,
    models.chosen = "all",
    em.by = 'all',
    em.algo = 'EMwmean',
    metric.select = 'ROC',
    metric.select.dataset = 'validation',
    metric.eval = c('TSS', 'ROC'),
    var.import = 3,
    EMci.alpha = 0.05,
    EMwmean.decay = "proportional",
  )
  
  ##Do projections
  sp_proj <-
    BIOMOD_Projection(
      bm.mod = myBiomodModelOut,
      proj.name = "predictors",
      new.env = get("bioclim_data", envir = globalenv()),
      models.chosen = 'all',
      metric.binary = 'ROC',
      build.clamping.mask = TRUE,
      do.stack = FALSE,
      output.format = ".tif"
    )
  
  ## ensemble model projections
  sp_ens_proj <- 
    BIOMOD_EnsembleForecasting(
      bm.em = myBiomodEM,
      bm.proj = sp_proj,
      models.chosen = 'all',
      metric.binary = 'ROC',
      compress = TRUE,
      do.stack = FALSE,
      output.format = ".tif"
    )
  return(paste0(sp," modelling completed !"))
}    

## launch the species modelling wrapper over species list ----
if(require(snowfall)){ ## parallel computation
  ## start the cluster
  sfInit(parallel = TRUE, cpus = 4) ## here we only require 4 cpus
  sfExportAll()
  sfLibrary(biomod2)
  ## launch our wrapper in parallel
  sf_out <- sfLapply(spp_to_model, biomod2_wrapper)
  ## stop the cluster
  sfStop()
} else { ## sequencial computation
  for (sp in spp_to_model){
    biomod2_wrapper(sp)
  }
  ## or with a lapply function in sequential model
  ## all_species_bm <- lapply(spp_to_model, biomod2_wrapper)
}      

load ("/Users/peggymutheu/Downloads/Ruppell's vulture/Gyps.rueppellii/Gyps.rueppellii.AllModels.ensemble.models.out")
get_evaluations(Gyps.rueppellii.AllModels.ensemble.models.out)

get_variables_importance(Gyps.rueppellii.AllModels.ensemble.models.out)
df <- print(get_variables_importance(Gyps.rueppellii.AllModels.ensemble.models.out))
df <- as.data.frame(df)
write.csv(df, "~gyps.rueppellii.csv")
Gyps_var <- read.csv("~gyps.rueppellii.csv", stringsAsFactors = FALSE)

my_sum <- Gyps_var %>%
  group_by(expl.var) %>%
  summarise( 
    n=n(),
    mean=mean(var.imp),
    sd=sd(var.imp)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

ggplot(my_sum) +
  geom_bar(aes(x=expl.var, y=mean), stat="identity",fill="darkgray", alpha=0.5) +
  geom_errorbar(aes(x=expl.var, ymin=mean-se, ymax=mean+se), width=0.5, colour="orange", alpha=0.6, size=0.6) +
  ggtitle("Ruppell's vulture")+ xlab("Variables ") + ylab("Relative Variable Importance")+ ylim(0, 0.4)+
  theme_classic()+
  coord_flip()+
  theme(plot.title = 
          element_text(hjust = 0.5))

get_built_models(Gyps.rueppellii.AllModels.models.out, full.name = NULL, PA = NULL, run = NULL, algo = NULL)
##Response curves for Ruppell's vulture
load("/Users/peggymutheu/Downloads/Ruppell's vulture/Gyps.rueppellii/Gyps.rueppellii.AllModels.ensemble.models.out")
get_built_models(Gyps.rueppellii.AllModels.ensemble.models.out)
mods <- get_built_models(Gyps.rueppellii.AllModels.ensemble.models.out, 
                         full.name = NULL,
                         merged.by.algo = 'mergedAlgo',
                         merged.by.run = 'mergedRun',
                         merged.by.PA = 'mergedData',
                         )
bm_PlotResponseCurves(bm.out = Gyps.rueppellii.AllModels.ensemble.models.out, 
                      models.chosen = mods,
                      fixed.var = 'median')
