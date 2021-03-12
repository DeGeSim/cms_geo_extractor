#pragma once

Cell *GeoExtractor::getCellptr(DetId &iterId)
{
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = getCellHashKeys(iterId);
  // Setup the detector
  Det &detector = detcol.detectors[detectorid];

  //Setup the subdetector
  Subdet &subdet = detector.subdetectors[subdetid];

  //Setup the layer
  Layer &layer = subdet.layers[layerid];

  Cell *cellptr;
  // setup the tiles or wafer depending on the detector
  if (isSiliconDet(detectorid))
  {
    //Setup the Wafer
    Wafer &wafer = layer.wafers[waferortileid];
    cellptr = &wafer.cells[cellid];
  }
  else
  {
    cellptr = &layer.tiles[waferortileid];
  }
  return cellptr;
}

CellHash GeoExtractor::getCellHashKeys(DetId &iterId)
{
  // Setup the detector
  unsigned int detectorid = iterId.det();
  Det &detector = detcol.detectors[detectorid];

  //Setup the subdetector
  unsigned int subdetid = iterId.subdetId();
  Subdet &subdet = detector.subdetectors[subdetid];

  //Setup the layer
  unsigned int layerid = recHitTools.getLayer(iterId);
  Layer &layer = subdet.layers[layerid];

  std::pair<int, int> waferortileid;
  std::pair<int, int> cellid;
  // setup the tiles or wafer depending on the detector
  if (isSiliconDet(detectorid))
  {
    //Setup the Wafer
    waferortileid = recHitTools.getWafer(iterId);
    cellid = recHitTools.getCell(iterId);
  }
  else
  {
    HGCScintillatorDetId scid = HGCScintillatorDetId(iterId);
    waferortileid = scid.ietaphi();
    cellid = std::make_pair(0, 0);
  }
  CellHash retuple = std::make_tuple(detectorid, subdetid, layerid, waferortileid, cellid);
  return retuple;
}



// this method is needed even though we can cout << cell
// because we need to print out cells that dont exitst.
std::string GeoExtractor::printcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> waferid, std::pair<int, int> cellid)
{
  std::ostringstream stringStream;
  stringStream << " Det " << detectorid;
  stringStream << " Subdet " << subdetid;
  stringStream << " Layer " << layerid;
  stringStream << " Wafer (" << waferid.first << "," << waferid.second << ")"<<"\n";
  stringStream << " Cell (" << cellid.first << "," << cellid.second << ")"<<"\n";
  std::string copyOfStr = stringStream.str();

  return copyOfStr;
}
