rm(list=ls())

setwd("~/Dropbox/SI_data_entry/assigning species")

sp.ids <- list(

  Vanessa_cardui = list(genus = "Vanessa", species = "cardui",
    subspecies = "", sex = "", author = "Linnaeus", family =
    "Nymphalidae", temp.id = c("JC_071612_2", "MM_080212_324",
    "PL_080912_353", "PL_081212_281", "PL_081112_283",
    "PL_081212_300", "PL_081112_266", "PL_081112_261",
    "PL_081112_262", "PL_081112_282", "PL_081112_284",
    "PL_081112_287", "PL_081212_299", "PL_081112_284B",
    "PL_081112_283B", "PL_081212_284", "PL_081012_5", "PL_080912_395",
    "PL_080912_397", "PL_081012_3", "MM_073112_125", "MM_073112_124",
    "MM_080212_224", "MM_073112_123", "MM_073112_130",
    "MM_080212_325", "MM_080412_250", "MM_080112_240",
    "PL_081312_291", "CH_082012_184", "PL_081312_280",
    "CH_082012_183", "PL_081312_228", "PL_080912_371",
    "PL_080912_410", "PL_080912_378", "PL_080912_396",
    "CH_082012_194", "PL_080912_351", "CH_082012_188",
    "MM_080312_211", "MM_080212_333", "MM_080212_332",
    "MM_080212_338", "MM_080412_246", "PL_081012_318",
    "PL_081312_290", "PL_081312_239", "PL_081012_338", "CH_081812_21",
    "PL_081312_294", "PL_081012_341", "CH_082012_187",
    "PL_081312_292", "PL_081012_4", "PL_080912_398", "PL_081312_234",
    "PL_081312_235", "PL_081012_323", "PL_081012_321",
    "PL_081012_322", "PL_080912_427", "PL_081312_227",
    "PL_081312_286", "PL_081312_284", "PL_081012_342",
    "MM_080112_126", "PL_081012_329", "MM_080112_249",
    "PL_081012_330", "PL_081012_326", "PL_081012_328",
    "PL_081012_327", "PL_081112_260", "PL_081312_232",
    "PL_081212_305", "JC_071612_6", "PL_081212_282", "PL_081212_283",
    "PL_081112_258", "JC_071412_23", "PL_081312_233", "PL_081112_259",
    "PL_081212_310", "CH_082012_193", "CH_082112_188",
    "PL_081112_257B", "CH_082112_203", "CH_081912_18",
    "CH_082112_198", "PL_081212_303B", "PL_081212_303C",
    "PL_081312_281", "PL_081312_289", "PL_081312_293")),

  Battus_philenor = list(genus = "Battus", species = "philenor",
    subspecies = "", sex = "", author = "Linnaeus", family =
    "Papilionidae", temp.id = c("PL_081212_307", "PL_081212_306",
    "PL_081212_278", "PL_081212_274", "PL_081312_223",
    "PL_081312_253", "PL_081312_224", "PL_081312_252",
    "PL_081312_222", "PL_081212_275", "PL_081112_252",
    "PL_081212_303", "PL_081212_304", "PL_081212_277")),

  Limenitis_bredowii = list(genus = "Limenitis", species = "bredowii",
    subspecies = "", sex = "", author = "Butler", family =
    "Nymphalidae", temp.id = c("CH_081812_20", "CH_082012_192",
    "CH_082012_163", "CH_082012_164", "CH_082112_184",
    "CH_082012_195")),

  Colias_eurytheme = list(genus = "Colias", species = "eurytheme",
    subspecies = "", sex = "", author = "Boisduval", family =
    "Pieridae", temp.id = c("SC_072712_265", "PL_081112_279",
    "PL_081012_325", "PL_080912_405", "PL_080912_383",
    "PL_081112_273", "PL_081112_267", "PL_081112_267B",
    "PL_081112_268", "PL_081212_288", "PL_081212_285",
    "PL_081112_288", "PL_081212_293", "PL_081212_292",
    "PL_081312_229", "CH_082112_189", "PL_080912_405",
    "PL_080912_384", "PL_080912_383", "PL_080912_391",
    "PL_080912_388", "PL_080912_381", "PL_080912_401",
    "PL_080912_400", "PL_080912_380", "PL_080912_411", "JC_071712_14",
    "PL_081212_294", "PL_081312_231", "PL_081112_278",
    "PL_081112_278B", "PL_081112_272", "PL_081112_275",
    "JC_071212_127", "PL_081112_277", "PL_081112_274",
    "PL_081212_291", "PL_081012_334", "CH_082112_199",
    "PL_081012_335", "CH_082112_186", "CH_082012_189",
    "CH_082112_187", "SC_072712_266", "SC_072712_268",
    "PL_081012_336", "PL_081012_325", "PL_081312_276",
    "MM_080212_335", "SC_072712_266", "SC_072712_263",
    "SC_072712_264", "SC_072612_283", "MM_080112_241",
    "MM_080212_331", "JC_071412_28", "JC_071612_9", "JC_071612_12",
    "JC_071212_129", "JC_071712_11", "JC_071612_152", "JC_071712_13",
    "JC_071612_10", "JC_071512_17", "JC_071212_126", "JC_071612_8",
    "JC_071512_9", "JC_071512_3", "JC_071612_3", "JC_071612_4",
    "JC_071712_6", "SC_072012_323", "JC_071512_16", "SC_072012_325",
    "PL_081012_333", "JC_071512_2", "JC_071512_18", "JC_071512_4",
    "JC_071612_11", "JC_071712_8", "JC_071612_5", "JC_071512_5",
    "SC_072012_327", "JC_071512_6", "SC_072712_267", "SC_072712_267B",
    "PL_081112_279")),

  Euptoierta_claudia = list(genus = "Euptoierta", species = "claudia",
    subspecies = "", sex = "", author = "Cramer", family =
    "Nymphalidae", temp.id = c("CH_082012_185", "CH_082112_185",
    "PL_081312_278", "PL_081312_236", "PL_081012_320",
    "JC_071212_128", "PL_081212_279", "CH_082012_182",
    "PL_081312_287", "PL_081312_288", "PL_081312_238",
    "MM_080412_245", "PL_081212_290")),

  Vanessa_carye_annabella = list(genus = "Vanessa", species = "carye",
    subspecies = "annabella", sex = "", author = "Field", family =
    "Nymphalidae", temp.id = c("PL_081212_309", "CH_082012_197",
    "PL_081112_289", "JC_071512_15", "MM_073112_128")),

  Vanessa_atlanta = list(genus = "Vanessa", species = "atlanta",
    subspecies = "", sex = "", author = "Linnaeus", family =
    "Nymphalidae", temp.id = c("CH_082012_167")),

  Vanessa_virginiensis = list(genus = "Vanessa", species =
    "virginiensis", subspecies = "", sex = "", author = "Drury",
    family = "Nymphalidae", temp.id =c("JC_071412_24",
    "CH_082112_201", "CH_082012_198")),

  Papilio_polyxenes = list(genus = "Papilio", species = "polyxenes",
    subspecies = "", sex = "", author = "Fabricius", family =
    "Papilionidae", temp.id = c("PL_080912_352", "PL_081212_280",
    "PL_081112_271", "PL_081212_276", "CH_082212_88")),

  Pontia_protodice = list(genus = "Pontia", species = "protodice",
    subspecies = "", sex = "", author = "Boisduval & Leconte", family
    = "Pieridae", temp.id = c("JC_071612_13", "JC_071512_10",
    "JC_071412_27", "JC_071512_7", "JC_071512_8", "JC_071212_130",
    "JC_071412_26", "JC_071512_11", "JC_071612_153",
    "CH_082112_202")),

  Hesperia_pahaska = list(genus = "Hesperia", species = "pahaska",
    subspecies = "", sex = "", author = "Leussler", family =
    "Hesperiidae", temp.id = c("PL_081012_314", "PL_081012_311",
    "PL_080912_363", "PL_080912_426", "PL_081012_312",
    "PL_081012_309", "PL_081112_111", "PL_081012_302",
    "PL_080912_365", "PL_081112_131", "PL_081012_304",
    "PL_080912_425", "PL_080912_423", "PL_081212_230",
    "PL_081112_133", "PL_080912_366", "PL_080912_420",
    "PL_080912_393", "PL_081112_114", "PL_081112_132",
    "PL_080912_370", "PL_081112_115", "PL_080912_367",
    "PL_081012_306", "PL_080912_368", "PL_080912_360",
    "PL_081112_257", "PL_081212_232", "PL_081212_250",
    "PL_081112_253B", "PL_081112_254", "PL_081112_285",
    "PL_081212_238", "PL_081212_237", "PL_081312_271",
    "PL_081312_270", "PL_081312_245", "PL_081312_247",
    "PL_081312_262", "PL_081312_263")),

  Hesperia_comma = list(genus = "Hesperia", species = "comma",
    subspecies = "", sex = "", author = "Linnaeus", family =
    "Hesperiidae", temp.id = c("PL_081312_261", "PL_081312_269",
    "PL_081212_251", "PL_081212_231", "PL_081112_256",
    "PL_080912_422", "PL_080912_416", "PL_080912_364",
    "PL_080912_375", "PL_081212_229", "PL_080912_371B",
    "PL_080912_362", "PL_080912_424", "PL_080912_386",
    "PL_080912_369", "PL_080912_374", "PL_080912_357",
    "PL_081112_270", "PL_080912_417", "PL_080912_361",
    "PL_080912_421", "PL_080912_418", "PL_080912_385",
    "PL_081212_245", "PL_081212_247", "PL_081212_249",
    "PL_081112_112", "PL_081112_130", "PL_081112_113",
    "PL_081012_313", "PL_081112_108", "PL_081012_310",
    "PL_081012_297", "PL_081012_300", "PL_081012_294",
    "PL_081012_295", "PL_081012_296", "PL_081112_110",
    "PL_081112_109", "PL_081312_250", "PL_081012_301",
    "PL_081012_298", "PL_081012_299")),

  Nathalis_iole = list(genus = "Nathalis", species = "iole",
    subspecies = "", sex = "", author = "Boisduval", family =
    "Pieridae", temp.id = c("MM_080212_320", "MM_080212_323",
    "SC_072812_107", "JC_071712_187", "JC_071612_7", "JC_071512_14")),

  Eurema_nicippe = list(genus = "Eurema", species = "nicippe",
    subspecies = "", sex = "", author = "Cramer", family = "Pieridae",
    temp.id = c("PL_081212_295", "PL_081312_275", "JC_071712_10",
    "JC_071212_131")),

  Polygonia_gracilis = list(genus = "Polygonia", species = "gracilis",
    subspecies = "", sex = "", author = "Grote & Robinson", family =
    "Nymphalidae", temp.id =c("MM_073112_122", "MM_080212_337",
    "MM_080212_336", "MM_080112_250", "CH_082012_165",
    "CH_082012_166", "MM_080212_334", "MM_080412_247",
    "MM_080212_328", "MM_080212_327", "SC_072912_123")),

  Plebejus_acmon = list(genus = "Plebejus", species = "acmon",
    subspecies = "", sex = "", author = "Westwood", family =
    "Lycaenidae", temp.id =c("SC_072712_270", "SC_072612_285",
    "SC_072012_329", "SC_072012_328")),

  Ministrymon_leda = list(genus = "Ministrymon", species = "leda",
    subspecies = "", sex = "", author = "Edwards", family =
    "Lycaenidae", temp.id =c("SC_072012_330")),

  Strymon_melinus = list(genus = "Strymon", species = "melinus",
    subspecies = "", sex = "", author = "Hubner", family =
    "Lycaenidae", temp.id =c("MM_073112_126", "JC_071712_7",
    "PL_081212_302", "JC_071712_186")),

  Leptotes_marina = list(genus = "Leptotes", species = "marina",
    subspecies = "", sex = "", author = "Reakirt", family =
    "Lycaenidae", temp.id =c("PL_081112_253", "CH_082012_201",
    "MM_080312_210", "PL_080912_415", "CH_082112_212")),

  Plebejus_glandon = list(genus = "Plebejus", species = "glandon",
    subspecies = "", sex = "", author = "de Prunner", family =
    "Lycaenidae", temp.id =c("JC_071712_4")),

  Celastrina_argiolus_ladon = list(genus = "Celastrina", species =
    "argiolus", subspecies = "ladon", sex = "", author = "Cramer",
    family = "Lycaenidae", temp.id = c("PL_081212_289",
    "MM_080312_213")),

  Phyciodes_tharos= list(genus = "Phyciodes", species = "tharos",
    subspecies = "", sex = "", author = "Drury", family =
    "Nymphalidae", temp.id =c("PL_081312_243", "PL_081012_339",
    "CH_082012_186", "PL_080912_387", "PL_080912_414")),

  Pyrgus_commonis = list(genus = "Pyrgus", species = "commonis",
    subspecies = "", sex = "", author = "Grote", family =
    "Hesperiidae", temp.id =c("PL_080912_419", "PL_081012_303",
    "CH_082012_196", "PL_080912_359", "PL_081212_236",
    "PL_081212_239", "JC_071512_119", "SC_072612_286",
    "JC_071512_173")),

  Ochlodes_snowi = list(genus = "Ochlodes", species = "snowi",
    subspecies = "", sex = "", author = "Edwards", family =
    "Hesperiidae", temp.id = c("PL_081212_248", "PL_081312_249")),

  Autochton_cellus = list(genus = "Autochton", species = "cellus",
    subspecies = "", sex = "", author = "Boisduval & Leconte", family
    = "Hesperiidae", temp.id =c("PL_081312_237")),

  Poanes_zabulon_taxiles = list(genus = "Poanes", species = "zabulon",
    subspecies = "taxiles", sex = "", author = "Edwards", family =
    "Hesperiidae", temp.id =c("SC_072012_324")),

  Epargyreus_clarus = list(genus = "Epargyreus", species = "clarus",
    subspecies = "", sex = "", author = "Cramer", family =
    "Hesperiidae", temp.id =c("PL_081012_324", "PL_081112_106",
    "PL_081112_281", "PL_081112_280", "PL_081212_297",
    "PL_081212_296", "PL_081212_298", "PL_080912_372",
    "PL_080912_406", "PL_081112_254", "PL_081312_279",
    "PL_081312_240", "PL_081112_256", "PL_081112_255")),

  Achalarus_casica = list(genus = "Achalarus", species = "casica",
    subspecies = "", sex = "", author = "Herrich-Schaffer", family =
    "Hesperiidae", temp.id =c("PL_081312_285", "PL_081312_244",
    "PL_081312_277","PL_081312_285")),

  Papilio_multicaudata = list(genus = "Papilio", species =
    "multicaudata", subspecies = "", sex = "", author = "Kirby",
    family = "Papilionidae", temp.id =c("PL_081112_251")),

  Neophasia_menapia = list(genus = "Neophasia", species = "menapia",
    subspecies = "", sex = "", author = "Felder & Felder", family =
    "Pieridae", temp.id =c("SC_072012_326")),

  Libytheana_carinenta = list(genus = "Libytheana", species =
    "carinenta", subspecies = "", sex = "", author = "Cramer", family
    = "Nymphalidae", temp.id =c("PL_081312_274")),

  Speyeria_mormonia = list(genus = "Speyeria", species = "mormonia",
    subspecies = "", sex = "", author = "Boisduval", family =
    "Nymphalidae", temp.id =c("JC_071712_2")),

  Cercyonis_oetus = list(genus = "Cercyonis", species = "oetus",
    subspecies = "", sex = "", author = "Boisduval", family =
    "Nymphalidae", temp.id =c("JC_071512_19")),

  Coenonympha_tullia = list(genus = "Coenonympha", species = "tullia",
    subspecies = "", sex = "", author = "Muller", family =
    "Nymphalidae", temp.id =c("JC_071612_1", "JC_071712_9")),

  Speyeria_atlantis = list(genus = "Speyeria", species = "atlantis",
    subspecies = "", sex = "", author = "Edwards", family =
    "Nymphalidae", temp.id = c("SC_072612_284", "PL_080912_349",
    "PL_080912_350", "PL_081012_331", "PL_081312_225",
    "PL_081012_343", "PL_080912_390", "PL_080912_354",
    "PL_080912_389", "PL_080912_402", "PL_080912_404",
    "PL_080912_376", "PL_081112_276", "PL_081312_272",
    "PL_081312_282", "PL_081312_283", "PL_081312_226", "JC_071712_1",
    "SC_072612_282", "PL_081012_332", "PL_081012_340",
    "PL_081212_301", "PL_081212_308", "PL_081312_230",
    "PL_081112_265", "PL_081012_315", "PL_080912_377",
    "PL_080912_382", "PL_081012_317", "PL_081012_319",
    "PL_081012_316", "PL_081012_337", "PL_080912_408",
    "PL_080912_409", "PL_080912_407", "PL_080912_399",
    "PL_080912_403")),


  Hemiargus_isola = list(genus = "Hemiargus", species = "isola",
    subspecies = "", sex = "", author = "Reakirt", family = "Lycaenidae",
    temp.id = c("SC_072812_106", "SC_072912_122", "PL_081312_241",
    "MM_080212_321", "MM_080412_248", "PL_081312_242",
    "MM_080212_339", "JC_071712_3", "CH_082212_90", "MM_080212_330",
    "MM_080212_329", "CH_082212_91", "CH_082012_191", "SC_072912_120",
    "CH_082012_199", "CH_082212_92", "SC_072912_121", "SC_072312_298",
    "SC_072812_108", "CH_082012_190", "MM_080112_248",
    "MM_080112_242", "MM_080212_326", "MM_080312_212",
    "MM_080212_340", "PL_080912_394", "CH_082112_209",
    "CH_082112_206", "CH_082112_208", "PL_080912_413",
    "CH_082112_210", "JC_071712_185", "PL_080912_412",
    "CH_082112_204", "PL_081112_263", "CH_082212_89", "CH_082112_211",
    "CH_082112_213", "PL_081112_264", "CH_082112_205")),

  Garisma_garita = list(genus = "Garisma", species = "garita",
    subspecies = "", sex = "", author = "Reakirt", family =
    "Hesperiidae", temp.id = c("MM_080112_245", "MM_080112_246",
    "MM_080112_247", "MM_080312_214", "MM_080312_215",
    "MM_080312_216", "MM_080312_217", "MM_080412_252",
    "MM_073112_127", "JC_071712_5", "JC_071712_181", "JC_071612_58",
    "JC_071212_124", "JC_071712_22", "JC_071212_125", "JC_071712_29",
    "JC_071512_174", "JC_071612_139", "JC_071612_142")),

  Hemiargus_ceraunus = list(genus = "Hemiargus", species = "ceraunus",
    subspecies = "", sex = "", author = "Fabricius", family = "Lycaenidae",
    temp.id = c("CH_082112_207")),

  Polites_draco = list(genus = "Polites", species = "draco",
    subspecies = "", sex = "", author = "Edwards", family =
    "Hesperiidae", temp.id = c("JC_071612_143")),

  Texola_elada = list(genus = "Texola", species = "elada", subspecies
    = "", sex = "", author = "Hewitson", family = "Nymphalidae",
    temp.id = c("PL_081112_269")),

  Atrytonopsis_morpho1 = list(genus = "Atrytonopsis", species = "sp.",
    subspecies = "a", sex = "", author = "", family = "Hesperiidae",
    temp.id = c("JC_071412_29")),

  Poanes_zabulon_taxiles = list(genus = "Poanes", species = "zabulon",
    subspecies = "taxiles", sex = "", author = "Edwards", family = "Hesperiidae",
    temp.id = c("CH_082212_39")),

  Thorybes_pylades = list(genus = "Thorybes", species = "pylades",
    subspecies = "", sex = "", author = "Scudder", family = "Hesperiidae", temp.id =
    c("PL_081112_105", "PL_080912_355", "PL_080912_392",
    "PL_081112_104", "PL_081212_235")),

  Erynnis_morpho1 = list(genus = "Erynnis", species = "sp.",
    subspecies = "a", sex = "", author = "", family = "Hesperiidae",
    temp.id = c("PL_081312_244")),

  Erynnis_morpho2 = list(genus = "Erynnis", species = "sp.",
    subspecies = "b", sex = "", author = "", family = "Hesperiidae",
    temp.id = c("PL_081312_277"))

  )


spec <- read.csv(file="Specimens06122013.csv", as.is=TRUE)

## use a for loop to add in species, genus, and sex
for(i in 1:length(sp.ids)) {
  ind <- spec$temp.id %in% sp.ids[[i]]$temp.id
  spec$species[ind] <- sp.ids[[i]]$species
  spec$subspecies[ind] <- sp.ids[[i]]$subspecies
  spec$genus[ind] <- sp.ids[[i]]$genus
  spec$sex[ind] <- sp.ids[[i]]$sex
  spec$author[ind] <- sp.ids[[i]]$author
  spec$determiner[ind] <- "L. Ponisio"
  spec$dateDetermined[ind] <- "2013"
  spec$order[ind] <- "Lepidoptera"
}

write.csv(spec, file="Specimens06122013.csv", row.names=FALSE)
