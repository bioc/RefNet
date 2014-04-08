# go.R
#------------------------------------------------------------------------------------------------------------------------
options(stringsAsFactors=FALSE)
library(RUnit)
source("../standardColumns.R")
#------------------------------------------------------------------------------------------------------------------------
run = function (levels)
{
  if (0 %in% levels) {
    tbl <<- read.table("hypoxiaSignaling-2006-v2.tsv", header=TRUE, sep="\t", as.is=TRUE)
    } # 0

  if (1 %in% levels) {
    bidirectional <- rep(FALSE, nrow(tbl))
    provider <- rep("pouyssegur-2006", nrow(tbl))
    tbl <<- cbind(tbl, bidirectional=bidirectional, provider=provider)
    checkTrue(length(setdiff(standardColumns, colnames(tbl))) == 0)
    tbl <<- tbl[, standardColumns]  # put them in the right order
    } # 1

  if (2 %in% levels) {
    write.table(tbl, sep="\t", row.names=FALSE, quote=FALSE, file="hypoxiaSignaling-2006.tsv")
    } # 2

  if (3 %in% levels) {
    } # 3

  if (4 %in% levels) {
    } # 4

  if (5 %in% levels) {
    } # 5

  if (6 %in% levels) {
    } # 6

  if (7 %in% levels) {
    } # 7

  if (8 %in% levels) {
    } # 8

  if (9 %in% levels) {
    } # 9

  if (10 %in% levels) {
    } # 10

  if (11 %in% levels) {
    } # 11

  if (12 %in% levels) {
    } # 12

  if (13 %in% levels) {
    } # 13

  if (14 %in% levels) {
    } # 14

  if (15 %in% levels) {
    } # 15

  if (16 %in% levels) {
    } # 16

  if (17 %in% levels) {
    } # 17

  if (18 %in% levels) {
    } # 18

  if (19 %in% levels) {
    } # 19

  if (20 %in% levels) {
    } # 20


} # run
#------------------------------------------------------------------------------------------------------------------------
