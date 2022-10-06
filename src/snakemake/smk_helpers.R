library(yaml)

setClass("smk", representation(
  config = "list",
  wildcards = "list",
  input = "list"
))
