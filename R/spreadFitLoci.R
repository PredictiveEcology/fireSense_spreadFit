#' Start cells and target sizes of the fires the spread model is fitted to
#'
#' The fires are those with buffers in `fireBufferedListDT` (`fireSense_dataPrepFit` leaves out fires with too
#' little burnable area). Each fire's target size is its burned area on burnable pixels (column
#' `fireSenseUtils::burnableHaTxt`, from `fireSense_dataPrepFit`), not its polygon area: the fit's landscape has
#' burnable pixels only, so a simulated fire cannot reach a polygon's non-burnable part.
#'
#' @param ras the raster of the fit's pixel IDs (`sim$rasterToMatch`).
#' @param pts list named `year<year>` of the fires' start points (`sim$spreadFirePoints`).
#' @param fireBufferedListDT list named `year<year>` of `data.table`s with the fires' `ids`.
#' @return as `fireSenseUtils::makeLociList()`.
spreadFitLoci <- function(ras, pts, fireBufferedListDT) {
  sizeCol <- fireSenseUtils::burnableHaTxt
  if (!sizeCol %in% names(pts[[1]]))
    stop("spreadFirePoints has no '", sizeCol, "' column (each fire's area on burnable pixels): ",
         "update fireSense_dataPrepFit, which makes it")
  idCol <- grep("ID", names(pts[[1]]), value = TRUE) |> setdiff("GID")
  yrs <- intersect(names(pts), names(fireBufferedListDT))
  pts <- Map(p = pts[yrs], fb = fireBufferedListDT[yrs], function(p, fb) p[as.data.frame(p)[[idCol]] %in% fb$ids, ])
  makeLociList(ras = ras, pts = pts, idsCol = idCol, sizeCol = sizeCol, yearPrefix = fireSenseUtils::yearTxt)
}
