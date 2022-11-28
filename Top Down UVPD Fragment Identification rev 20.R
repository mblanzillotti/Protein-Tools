#-------------------------------------------------------------------------------
#Fragment Matching by Isotopic Distributions in m/z Space
#-------------------------------------------------------------------------------

#===============================================================================
#This scrippt is meant to provide matching of fragment ions produced in UVPD
#spectra of proteins or peptides in the m/z domain based on their isotopic
#distributions. Raw files from Thermo Orbitrap instruments can be uploaded
#directly, alongside a primary sequence with support for custom modifications
#endoded in parentheses of the form: numeric - (+XX.XXX), formula - (CXHYOZ), or
#predefined = (Acetyl).

#Built in R 4.1.3 "One Push-Up" - 23 March 2022
#===============================================================================


#-------------------------------------------------------------------------------
#Prerequisites
#-------------------------------------------------------------------------------

#Package Requirements-----------------------------------------------------------
require(tidyverse)
require(readxl)
require(plotly)
require(scales)
require(grid)
require(ggthemes)
require(IsoSpecR)
require(rawrr)

#Elemental Compositions---------------------------------------------------------
element_masses <- 
  tribble(
    ~element, ~mono_mass,
    "C",      12.00000,
    "H",      1.007825035,
    "N",      14.003074,
    "O",      15.99491603,
    "P",      30.973762,
    "S",      31.9720707,
    "Na",     22.9897677,
    "Zn",     63.929145,
    "Fe",     53.939612
  )
## Element monoisotopic masses.  If an element is not defined here, and in the following
## tribble, it will not be considered during theoretical mass calculation.

symbol_composition <- 
  tribble(
    ~symbol,    ~C, ~H, ~N, ~O, ~P, ~S, ~Fe,
    "A",        3,  5,  1,  1,  0,  0,  0,  
    "C",        3,  5,  1,  1,  0,  1,  0,
    "D",        4,  5,  1,  3,  0,  0,  0,
    'E',        5,  7,  1,  3,  0,  0,  0,  
    "F",        9,  9,  1,  1,  0,  0,  0,  
    "G",        2,  3,  1,  1,  0,  0,  0,  
    "H",        6,  7,  3,  1,  0,  0,  0,  
    "I",        6,  11, 1,  1,  0,  0,  0,  
    "K",        6,  12, 2,  1,  0,  0,  0,  
    "L",        6,  11, 1,  1,  0,  0,  0,  
    "M",        5,  9,  1,  1,  0,  1,  0,  
    "N",        4,  6,  2,  2,  0,  0,  0,  
    'P',        5,  7,  1,  1,  0,  0,  0,  
    "Q",        5,  8,  2,  2,  0,  0,  0,  
    "R",        6,  12, 4,  1,  0,  0,  0,  
    "S",        3,  5,  1,  2,  0,  0,  0,  
    "T",        4,  7,  1,  2,  0,  0,  0,  
    "V",        5,  9,  1,  1,  0,  0,  0,
    "W",        11, 10, 2,  1,  0,  0,  0,  
    "Y",        9,  9,  1,  2,  0,  0,  0,
    "Acetyl",   2,  2,  0,  1,  0,  0,  0,
    "Carbamyl", 1,  2,  1,  1,  0,  0,  0,
    "Methyl",   1,  2,  0,  0,  0,  0,  0,
    "Deoxid",   0,  0,  0, -1,  0,  0,  0,
    "Phospho",  0,  1,  0,  3,  1,  0,  0,
    "Succinyl", 4,  4,  0,  2,  0,  0,  0,
    "Nterm",    0,  1,  0,  0,  0,  0,  0,
    "Cterm",    0,  1,  0,  1,  0,  0,  0,
    "Heme",     34, 29, 4,  4,  0,  0,  1
  ) %>%
  pivot_longer(-symbol, names_to = "element", values_to = "amount")

## Predefined chemical formulae for defined elements. Includes amino acid residues and some PTM's 

#Parsing Functions--------------------------------------------------------------

interpret_formula <- function(encoding){
  elemental_composition <- tibble(
    element = str_extract_all(encoding, "[:alpha:]{1}[:lower:]?\\d+") %>% 
      as_vector() %>% 
      str_extract("[:alpha:]{1}[:lower:]?"),
    total = str_extract_all(encoding, "[:upper:]{1}[:lower:]?\\d+") %>% 
      as_vector() %>% 
      str_extract("\\d+") %>% 
      as.numeric()
  )
  return(elemental_composition)
}
## Identifies and interprets a chemical formula in parentheses during sequence parsing.

interpret_digits <- function(encoding){
  return(
    str_remove_all(encoding, "[\\(\\)]") %>% 
    str_extract("[\\-\\+]?\\d*\\.?\\d*") %>% 
      as.numeric()
  )
}
## Identifies and interprets a numerical modification in parentheses during sequence parsing.

interpret_symbol <- function(encoding){
  return(filter(symbol_composition,
                str_detect(symbol, paste("^", encoding, "$", sep = ""))
                )
  )
}
# Identifies a predefined chemical modification in parentheses during sequence parsing

composition_to_mass <- function(composition){
  return(
    composition %>%
      left_join(element_masses, by = "element") %>%
      summarise(sum(total*mono_mass)) %>%
      as_vector() %>%
      unname()
  )
}
## Converts a chemical formula of elements defined in "element_masses" from a composition to mass

composition_to_nominal <- function(composition){
  return(
    composition %>%
      left_join(element_masses, by = "element") %>%
      summarise(sum(total*mono_mass)) %>%
      as_vector() %>%
      unname() %>% 
      round(0)
  )
}
## Converts a chemical formula of elements defined in "element_masses" from a composition to an integer mass

check_termini <- function(input_seq){
  
  input_seq <- str_remove_all(input_seq, "\\s")
  
  if(str_detect(input_seq, "Nterm") && str_detect(input_seq, "Cterm")){
    
    return(
      input_seq
    )
    
  } else if(!str_detect(input_seq, "Nterm") && str_detect(input_seq, "Cterm")){
    
    return(
      paste("Nterm", input_seq, sep = "")
    )
    
  } else if(str_detect(input_seq, "Nterm") && !str_detect(input_seq, "Cterm")){
    
    return(
      paste(input_seq, "Cterm", sep = "")
    )
    
  } else {
    
    return(
      paste("Nterm", input_seq, "Cterm", sep = "")
    )
  }
  
}
## Inserts "Nterm" and "Cterm" on the N and C terminus of the sequence, if they are not there already, for clean parsing

parse_prot_input <- function(input_seq){

  input_seq <- check_termini(input_seq)

  parsed <- tibble(
    input = input_seq 
    %>% str_remove_all("\\s")
    %>% str_extract_all("[:upper:](term)?(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))?") %>% as_vector()
  ) %>% 
    mutate(
      residue = str_extract(input, "[:upper:](term)?"),
      modification = str_extract(input, "\\(\\+?\\-?[:alnum:]+[\\.]?[\\d]{0,}\\) ?") %>% str_remove_all("[\\(\\)]"),
      position = c(NA, 1:(length(input)-2), NA)
    ) %>% 
    select(-input)
  
  return(parsed)
}
## Parses the expected protein sequence input by residue, with any accompanying modifications

determine_composition <- function(parsed_input){
  
  residues <- left_join(parsed_input, symbol_composition, by = c("residue" = "symbol")) %>% 
    group_by(element) %>% 
    summarise(
      residue_total = sum(amount)
    )
  
  modifications <- filter(parsed_input, !is.na(modification)) %>% 
    left_join(symbol_composition, by = c("modification" = "symbol")) %>% 
    group_by(element) %>% 
    summarise(
     mod_total = sum(amount)
    ) %>% 
    filter(!is.na(element))
  
  formula_mods <- filter(parsed_input, !is.na(modification)) %>% 
    left_join(symbol_composition, by = c("modification" = "symbol")) %>% 
    filter(str_detect(modification, "[:alpha:]")) %>% 
    filter(str_detect(modification, "[:digit:]")) %>% 
    select(modification) %>% 
    as_vector()
    
  formula_comps <- tibble(
    element = vector("character"),
    total = vector("numeric"),
  )
  
  if(length(formula_mods) != 0){
    
    for(i in 1:length(formula_mods)){
      formula_comps <- add_row(formula_comps,
        interpret_formula(formula_mods[i])
      )
    }
    
  }
  
  
  formula_comps <- rename(formula_comps, formula_total = total)
    
    composite <- full_join(residues, modifications, by = "element") %>% 
      full_join(formula_comps, by = "element") %>% 
      mutate(
        formula_total = ifelse(is.na(formula_total), 0, formula_total),
        residue_total = ifelse(is.na(residue_total), 0, residue_total),
        mod_total = ifelse(is.na(mod_total), 0, mod_total),
        total = mod_total + residue_total + formula_total
      )
  
  return(composite %>% select(element, total))
}
## Yields a chemical formula of predefined elements from the parsed protein sequence


calc_monoiso <- function(parsed_input){
  
    composition_mass <- determine_composition(parsed_input) %>%
      composition_to_mass()
    
    modifications <- parsed_input %>%
      filter(!is.na(modification))
    
    if(nrow(modifications) == 0 || all(!str_detect(modifications$modification, "[\\-\\+]?\\d+\\.?\\d*"))){
      
      digit_mass <- 0
      
    } else if(any(str_detect(modifications$modification, "[\\-\\+]?\\d+\\.?\\d*"))){
      
      digit_mass <- filter(modifications, !str_detect(modification, "[:alpha:]")) %>%
        mutate(
          mass = vapply(modification, interpret_digits, numeric(1)) %>% unname()
        ) %>%
        select(mass) %>% 
        as_vector() %>% 
        sum()
      
    } else {
      
      digit_mass <- 0
      
    }
    
    composition_mass <- composition_mass + digit_mass
   
  
  return(composition_mass) 
}
## Calculates a monoisotopic mass from the set of predefined elements, adding any numerical inputs.

fragment_composition <- function(parsed_input, pos, ion_type){
  
  if(str_detect(ion_type, "^a$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Cterm")) %>% 
                              filter(position <= pos | str_detect(residue, "Nterm"))) %>% 
        mutate(total = ifelse(str_detect(element, "^C$"), total-1, total)) %>% 
        mutate(total = ifelse(str_detect(element, "^O$"), total-1, total))
    )
    
  } else if(str_detect(ion_type, "^a\\+1$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Cterm")) %>% 
                              filter(position <= pos | str_detect(residue, "Nterm"))) %>% 
        mutate(total = ifelse(str_detect(element, "^C$"), total-1, total)) %>%
        mutate(total = ifelse(str_detect(element, "^H$"), total+1, total)) %>%  
        mutate(total = ifelse(str_detect(element, "^O$"), total-1, total))
    )
    
  } else if(str_detect(ion_type, "^b$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Cterm")) %>% 
                              filter(position <= pos | str_detect(residue, "Nterm")))
    )
    
  } else if(str_detect(ion_type, "^b-H2O$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Cterm")) %>% 
                              filter(position <= pos | str_detect(residue, "Nterm"))) %>% 
        mutate(total = ifelse(str_detect(element, "^H$"), total-2, total)) %>% 
        mutate(total = ifelse(str_detect(element, "^O$"), total-1, total))
    )
    
  } else if(str_detect(ion_type, "^c$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Cterm")) %>% 
                              filter(position <= pos | str_detect(residue, "Nterm"))) %>% 
        mutate(total = ifelse(str_detect(element, "^H$"), total+3, total)) %>% 
        mutate(total = ifelse(str_detect(element, "^N$"), total+1, total))
    )
    
  } else if(str_detect(ion_type, "^x$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Nterm")) %>% 
                              filter(position > pos | str_detect(residue, "Cterm"))) %>% 
        mutate(total = ifelse(str_detect(element, "^C$"), total+1, total)) %>% 
        mutate(total = ifelse(str_detect(element, "^O$"), total+1, total))
    )
    
  } else if(str_detect(ion_type, "^x\\+1$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Nterm")) %>% 
                              filter(position > pos | str_detect(residue, "Cterm"))) %>% 
        mutate(total = ifelse(str_detect(element, "^C$"), total+1, total)) %>%
        mutate(total = ifelse(str_detect(element, "^H$"), total+1, total)) %>%
        mutate(total = ifelse(str_detect(element, "^O$"), total+1, total))
    )
    
  } else if(str_detect(ion_type, "^y$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Nterm")) %>% 
                              filter(position > pos | str_detect(residue, "Cterm"))) %>% 
        mutate(total = ifelse(str_detect(element, "^H$"), total+1, total))
    )
    
  } else if(str_detect(ion_type, "^y\\-1$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Nterm")) %>% 
                              filter(position > pos | str_detect(residue, "Cterm"))) 
    )
    
  } else if(str_detect(ion_type, "^y\\-2$") && str_detect(parsed_input$residue[pos+2], "P")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Nterm")) %>% 
                              filter(position > pos | str_detect(residue, "Cterm")))%>% 
        mutate(total = ifelse(str_detect(element, "^H$"), total-1, total))
    )
    
  } else if(str_detect(ion_type, "^z$")){
    
    return(
      determine_composition(parsed_input %>%
                              filter(!str_detect(residue, "Nterm")) %>% 
                              filter(position > pos | str_detect(residue, "Cterm"))) %>% 
        mutate(total = ifelse(str_detect(element, "^H$"), total-1, total)) %>% 
        mutate(total = ifelse(str_detect(element, "^N$"), total-1, total))
    )
    
  } else {
    
    return("Please input a valid ion type")
    
  }
  
}
## Generates a chemical formula for a single fragment ion type at a single position.

calc_frag_formula <- function(parsed_input, ion_types = c("a", "a+1", "b", "c", "x", "x+1", "y", "y-1", "y-2", "z")){
  
  parsed_input <- parsed_input %>% 
    mutate(
      digit_mods = map(modification, interpret_digits) %>% ifelse(is.na(.), 0, .),
      fwd_cumsum = cumsum(digit_mods),
      back_cumsum = rev(cumsum(rev(digit_mods)))
    )
  
  fragment_compositions <- tibble(
    ion_type = vector("character"),
    position = vector("integer"),
    digit_mass = vector("numeric"),
    element = vector("character"),
    total = vector("integer")
  )
  
  for(i in 1:(max(filter(parsed_input, !is.na(position))$position)-1)){
    
    for(j in seq_along(ion_types)){
      
      carrier <- fragment_composition(
        parsed_input,
        i,
        ion_types[j]
      )
      
      if(typeof(carrier) == "character"){
        
        next
        
      } else {
        
        fragment_compositions <- fragment_compositions %>% 
          add_row(
            ion_type = ion_types[j],
            position = i,
            digit_mass = ifelse(str_detect(ion_type, "[abc]"),
                                parsed_input$fwd_cumsum[i+1],
                                parsed_input$back_cumsum[i+1]
            ),
            element = carrier$element,
            total = carrier$total
            
          )
          
      }
      
    }
  }

  fragment_compositions <- fragment_compositions %>% 
    group_by(ion_type, position) %>% 
    nest(data = c(element, total)) %>%
    mutate(
      formula_mass = map(data, composition_to_mass) %>% as_vector(),
      calc_mass = formula_mass + digit_mass
    ) %>% 
    select(-digit_mass)

  return(fragment_compositions)
}
## Iterates "fragment_composition" over all positions and requested ion types to generate
## chemmical formulae for each fragment ion along the protein backbone.

composition_to_IsoSpec <- function(composition){
  IsoSpec_comp <- composition$total %>% as_vector()
  names(IsoSpec_comp) <- composition$element %>% as_vector()
  
  return(IsoSpec_comp)
} 
## Converts chemical formulae tibbles into a named vector format that is friendly with IsoSpec

#Isotope Generation Functions-------------------------------------------------

IsoSpec_nominal_dist <- function(IsoSpec_comp){
  
  nominal_dist <- IsoSpecR::IsoSpecify(IsoSpec_comp, 0.95) %>% 
    as_tibble() %>% 
    mutate(
      nominal_mass = round(mass, 0)
    ) %>% 
    group_by(nominal_mass) %>% 
    summarise(
      nominal_prob = sum(prob)
    ) %>% 
    select(nominal_mass, nominal_prob) %>% 
    slice_max(order_by = nominal_mass, n = 15) %>% 
    ungroup()
  
  return(nominal_dist)
    
}
## Generates a fine isotopic distribution using IsoSpec, and condenses fine values based on nominal mass

generate_iso_dist <- function(fragment_compositions){
  
  fragment_distributions <- fragment_compositions %>% 
    mutate(
      nominal_formula = round(formula_mass, 0),
      nominal_isos = map(data, composition_to_IsoSpec) %>%
        map(IsoSpec_nominal_dist)
    ) %>% 
    unnest(cols = c(nominal_isos)) %>%
    mutate(
      exact_mass = (nominal_mass - nominal_formula)*1.008665 + formula_mass + (calc_mass - formula_mass)
    ) %>% 
    select(ion_type, position, nominal_prob, exact_mass) %>% 
    nest(exact_isotopes = c(nominal_prob, exact_mass))
  
  return(fragment_distributions)
}
## Iterates "IsoSpec_nominal_dist" along fragment ion formulae to generate isotopic patterns for
## all fragments

charge_convolve <- function(fragment_distributions, min_charge = 1, max_charge, low_mz = 200, high_mz = 3000){
  
  fragment_mz <- fragment_distributions %>% 
    full_join(tibble(z = min_charge:max_charge), by = character(0)) %>% 
    unnest(cols = exact_isotopes) %>% 
    mutate(
      mz = (exact_mass + z*1.00727647)/abs(z)
    ) %>% 
    ungroup() %>% 
    group_by(ion_type, position, z) %>% 
    filter(mz > low_mz) %>% 
    filter(mz < high_mz) %>% 
    nest(iso_mz = c(exact_mass, mz, nominal_prob)) %>% 
    ungroup()
  
  return(fragment_mz)
}
## Convolves fragment ion distributions in m/z space with the requested number of charges
## and within the requested m/z range

ID_zcon_frags <- function(fragments, centroids, tolerance = 0.00001){
  
  matched_fragments <- tibble(
    ion_type = vector("character"),
    position = vector("integer"),
    charge = vector("integer"),
    exact_mass = vector("numeric"),
    mz = vector("numeric"),
    nominal_prob = vector("numeric"),
    expr_mz = vector("numeric"),
    intensity = vector("numeric"),
    mz_error = vector("numeric"),
    ppm_error = vector("numeric"),
    predicted_i = vector("numeric")
  )
  
  for(i in 1:nrow(fragments)){
    
    current_fragment <- fragments %>% 
      slice(i) %>% 
      unnest(iso_mz)
    
    isotope_matches <- tibble(
      ion_type = vector("character"),
      position = vector("integer"),
      charge = vector("integer"),
      exact_mass = vector("numeric"),
      mz = vector("numeric"),
      nominal_prob = vector("numeric"),
      expr_mz = vector("numeric"),
      intensity = vector("numeric"),
      series_no = vector("integer")
    )
    
    for(j in 1:nrow(current_fragment)){
      
      ID <- filter(centroids, near(centroids$expr_mz, current_fragment$mz[j], tol = current_fragment$mz[j]*tolerance))
      
      if(nrow(ID) > 0 && min(ID$intensity) > 10){
        
        ID <- mutate(ID,
                     mz_error = current_fragment$mz[j] - expr_mz,
              ) %>% 
              filter(abs(mz_error) == min(abs(mz_error)))
        
        isotope_matches <- add_row(isotope_matches,
          ion_type = current_fragment$ion_type[j],
          position = current_fragment$position[j],
          charge = current_fragment$z[j],
          exact_mass = current_fragment$exact_mass[j],
          mz = current_fragment$mz[j],
          nominal_prob = current_fragment$nominal_prob[j],
          expr_mz = ID$expr_mz,
          intensity = ID$intensity,
          series_no = j
        )
        
      }
      
    }
    
    #Diagnostic to evaluate ppm matching
    # if(sum(isotope_matches$nominal_prob) > 0.6){
    # 
    #   matched_fragments <- add_row(matched_fragments,
    #                                  select(isotope_matches, -series_no),
    #                                  predicted_i = sum(isotope_matches$intensity)*nominal_prob,
    #                                  mz_error = mz - expr_mz,
    #                                  ppm_error = mz_error / mz * 1e6
    #                        )
    # 
    # 
    # }
    
    if(nrow(isotope_matches) > 0){

      if(sum(isotope_matches$nominal_prob) > 0.60){

        isotope_fit_check <- mutate(isotope_matches,
          normalized_intensity = intensity / max(intensity),
          normalized_prob = nominal_prob / max(nominal_prob),
          frac_intensity = intensity / sum(intensity),
          prob_diff = abs(normalized_prob - normalized_intensity),
          frac_diff = nominal_prob - frac_intensity,
          mz_error = mz - expr_mz,
          ppm_error = mz_error / mz * 1e6,
          series_no = series_no,
          predicted_i = nominal_prob*sum(intensity),
          predicted_i_diff_pct = (predicted_i - intensity) / ((predicted_i + intensity)/2)
        )

        if(all(abs(isotope_fit_check$frac_diff) < 0.30) &&
           all(abs(diff(isotope_fit_check$series_no)) == 1) &&
           (sum(abs(isotope_fit_check$predicted_i_diff_pct) > 1) <= 3)
           ){

          matched_fragments <- add_row(matched_fragments,
              select(isotope_matches, -series_no),
              select(isotope_fit_check, predicted_i, mz_error, ppm_error)
          )

        }

      }

    }
    
  }
  
  
  return(matched_fragments)
}
## Identifies fragment ion isotopic distributions based on a requested ppm error
## and empirically defined tolerances.

#Graphical Functions----------------------------------------------------------

frag_colors <- c(
  "a" = "#008000",
  "a+1" = "#008000",
  "b" = "#0000FF",
  "c" = "#FF0000",
  "x" = "#008000",
  "x+1" = "#008000",
  "y" = "#0000FF",
  "y-1" = "#0000FF",
  "y-2" = "#0000FF",
  "z" = "#FF0000"
)

#-------------------------------------------------------------------------------
#Usage
#-------------------------------------------------------------------------------

date()#Start time

protein <- "Nterm(Acetyl)GDVEKGKKIFVQKC(-1H)AQC(Heme)HTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"
#Input protein sequence, in this case, ubiquitin

setwd("D:/UT Austin/Research/Top-Down/Publication Data/CYC")
#Set working directory, basically all the directories up to where you keep files

mass_list <- readSpectrum("07242022_CYC_13+_951_213nm_20ms_8e5_240k.raw", 1)
#Read in rawfile including metadata

picked_peaks <- tibble(
  expr_mz = mass_list[[1]]$centroid.mZ,
  intensity = mass_list[[1]]$centroid.intensity
)
#Put centroids from rawfile into a function-friendly tibble

z_convolved_frags <- parse_prot_input(protein) %>% 
  calc_frag_formula(ion_types = c("a", "a+1", "b", "c", "x", "x+1", "y", "y-1", "y-2", "z")) %>% 
  generate_iso_dist() %>% 
  charge_convolve(1, 13, low_mz = 200, high_mz = 2000) #%>% 
  # unnest(iso_mz) %>% 
  # filter(near(mz, 694, 1)) %>% 
  # filter(z == 4) %>% 
  # nest(iso_mz = c(exact_mass, mz, nominal_prob))
#Generate theoretical fragment ion distributions in m/z space based on input sequence

frag_IDs <- ID_zcon_frags(z_convolved_frags, picked_peaks, tolerance = 1.0e-5)
#Match fragment ions and store in a big ol tibble

date()#End time, if you're interested in that sort of thing

#-------------------------------------------------------------------------------
#Plotting Results
#-------------------------------------------------------------------------------
#I've left a helpful plotting functions to visualize results. Happy Hunting

#Fragment identifications should be tidy, which enables straightforward plotting
#and summarisation of fragmentation data.

#-------------------------------------------------------------------------------

plot_fragment_IDs <- function(rawfile, IDs){
  
  spectrum_profile <- tibble(
    mz = rawfile[[1]]$mZ,
    intensity = rawfile[[1]]$intensity
  )
  
  frag_colors <- c(
    "a" = "#008000",
    "a+1" = "#008000",
    "b" = "#0000FF",
    "c" = "#FF0000",
    "x" = "#008000",
    "x+1" = "#008000",
    "y" = "#0000FF",
    "y-1" = "#0000FF",
    "y-2" = "#0000FF",
    "z" = "#FF0000"
  )
  
  return(
    ggplot()+
      geom_line(
        data = spectrum_profile,
        aes(x = mz, y = intensity)        
      )+
      geom_col(
        data = IDs,
        aes(x = expr_mz, y = predicted_i, fill = ion_type)
      )+
      geom_point(
        data = IDs,
        aes(x = mz, y = predicted_i + 25, color = ion_type)
      )+
      theme_publication()+
      scale_color_manual(values = frag_colors)+
      scale_fill_manual(values = frag_colors)+
      ylab("Intensity")+
      xlab("m/z")
  )
}
#outputs an interactive, annotated spectrum
ggplotly(
  plot_fragment_IDs(mass_list, frag_IDs)+
    geom_col(data = picked_peaks, aes(x = expr_mz, y = intensity))
)

#-------------------------------------------------------------------------------

plot_fragment_histogram <- function(IDs){
  
  frag_colors <- c(
    "a" = "#008000",
    "a+1" = "#008000",
    "b" = "#0000FF",
    "c" = "#FF0000",
    "x" = "#008000",
    "x+1" = "#008000",
    "y" = "#0000FF",
    "y-1" = "#0000FF",
    "y-2" = "#0000FF",
    "z" = "#FF0000"
  )
  
  return(
    ggplot()+
      geom_col(
        data = IDs,
        aes(x = position, y = intensity, fill = ion_type), position = "stack"
      )+
      theme_publication()+
      scale_fill_manual(values = frag_colors)+
      ylab("Intensity")+
      xlab("Cleavage Position")
  )
}
#outputs a histogram where identified ion intensity is plotted against position
plot_fragment_histogram(frag_IDs)

#-------------------------------------------------------------------------------

plot_fragment_counts <- function(IDs){
  frag_colors <- c(
    "a" = "#008000",
    "a+1" = "#008000",
    "b" = "#0000FF",
    "c" = "#FF0000",
    "x" = "#008000",
    "x+1" = "#008000",
    "y" = "#0000FF",
    "y-1" = "#0000FF",
    "y-2" = "#0000FF",
    "z" = "#FF0000"
  )
  
  IDs <- mutate(IDs,
                  ion_family = str_extract(ion_type, "[abcxyz]")
                ) %>% 
    select(position, ion_family) %>% 
    distinct()
  
  return(
    ggplot(IDs)+
      geom_bar(aes(x = ion_family, fill = ion_family))+
      theme_publication()+
      scale_fill_manual(values = frag_colors)+
      xlab("Ion Type")+
      ylab("Number of Identified Fragments")+
      labs(fill = "Ion Type")
  )
  
}
#Summarizes number of unique fragment ions identified by ion type
plot_fragment_counts(frag_IDs)

#-------------------------------------------------------------------------------

plot_fragment_sum_i <- function(IDs){
  frag_colors <- c(
    "a" = "#008000",
    "a+1" = "#008000",
    "b" = "#0000FF",
    "c" = "#FF0000",
    "x" = "#008000",
    "x+1" = "#008000",
    "y" = "#0000FF",
    "y-1" = "#0000FF",
    "y-2" = "#0000FF",
    "z" = "#FF0000"
  )
  
  IDs <- mutate(IDs,
                ion_family = str_extract(ion_type, "[abcxyz]")
  ) %>% 
    select(ion_family, intensity, expr_mz, charge) %>% 
    distinct()
    
  
  return(
    ggplot(IDs)+
      geom_col(aes(x = ion_family, y = intensity, fill = ion_family))+
      theme_publication()+
      scale_fill_manual(values = frag_colors)+
      xlab("Ion Type")+
      ylab("Sum Intensity\nof Identified Fragments")+
      labs(fill = "Ion Type")
  )
  
}
#Summarizes total intensity of identified fragments by ion type
plot_fragment_sum_i(frag_IDs)

#-------------------------------------------------------------------------------

calc_sequence_cov <- function(protein, IDs){
  return(
    IDs$position %>%
      as_vector() %>% 
      unique() %>% 
      length() / (nrow(parse_prot_input(protein)) - 3)
  )
}
#Calculates sequence coverage.  Pretty straightforward
calc_sequence_cov(protein, frag_IDs)
