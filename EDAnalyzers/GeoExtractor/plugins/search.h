#pragma once

// wraps findNextCell and loops over the ids
void GeoExtractor::assignZneighbors(std::vector<DetId> &v_validHGCalIds)
{
  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    DetId iterId = v_validHGCalIds[i];
    // if (iterId.rawId() != 2426030080)
    // {
    //   continue;
    // }
    // LOG(DEBUG) << "skipping everything but 2426030080"<<"\n";
    if (i % 1000 == 0)
    {
      LOG(DEBUG) << "Assinging z neighbors " << i << "\n";
    }
    Cell *cellptr = getCellptr(iterId);
    validateId(cellptr->globalid);
    cellptr->next = findNextCell(cellptr->globalid);
    cellptr->previous = findPreviousCell(cellptr->previous);
  }
}

// This function decices in which layer in which detector(s) to search for the neigbor.
// The the act lifting is done by the findCellCloseToXYpos
DetId GeoExtractor::findNextCell(DetId cellId)
{
  LOG(DEBUG) << "Start findNextCell"
             << "\n";
  CellHash hash = getCellHashKeys(cellId);
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
  LOG(DEBUG) << "find cell for id" << cellId.rawId() << "\n";
  LOG(DEBUG) << getCellHashKeys(cellId) << "\n";
  //
  // HGCalEE = 8, layer 1-28
  // HGCalHSi = 9, layer 1-22
  // HGCalHSc = 10, layer 9-22
  // HGCalTrigger = 11 X

  // For the ee cal we can easily move forward
  if (detectorid == DetId::HGCalEE)
  {
    if (layerid < 28)
    {
      LOG(DEBUG) << "A"
                 << "\n";
      return findCellCloseToXYpos(cellId, detectorid, subdetid, layerid + 1).first;
    }
    // to the the next detector EE-> HSi
    else
    {
      LOG(DEBUG) << "B"
                 << "\n";
      // TODO check if there are nodes in layer 1 of the HGCalSc
      return findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, 1).first;
    }
  }
  if (detectorid == DetId::HGCalHSi || detectorid == DetId::HGCalHSc)
  {
    // for layer <8 all cells we can just search in the Si part
    if (layerid < 8 && detectorid == DetId::HGCalHSi)
    {
      LOG(DEBUG) << "C"
                 << "\n";
      return findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, layerid + 1).first;
    }
    // return 0 in the last layer
    else if (layerid == 22)
    {
      LOG(DEBUG) << "D"
                 << "\n";
      return DetId(0);
    }
    else
    {
      LOG(DEBUG) << "E"
                 << "\n";
      auto [sicanid, deltasi] = findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, layerid + 1);
      auto [sccanid, deltasc] = findCellCloseToXYpos(cellId, DetId::HGCalHSc, subdetid, layerid + 1);
      if (deltasi < deltasc)
      {
        return sicanid;
      }
      else
      {
        return sccanid;
      }
    }
  }
  // LOG(ERROR)  << hash << ""<<"\n"; //doenst work ?
  std::cout << hash;
  LOG(ERROR) << "Wont find neighbor. This part should never be reached."
             << "\n";
  exit(EXIT_FAILURE);
  return DetId(0);
}

DetId GeoExtractor::findPreviousCell(DetId cellId)
{
  LOG(DEBUG) << "Start findPreviousCell"
             << "\n";
  CellHash hash = getCellHashKeys(cellId);
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
  LOG(DEBUG) << "find cell for id" << cellId.rawId() << "\n";
  LOG(DEBUG) << hash << "\n";
  //
  // HGCalEE = 8, layer 1-28
  // HGCalHSi = 9, layer 1-22
  // HGCalHSc = 10, layer 9-22
  // HGCalTrigger = 11 X

  // For the ee cal we can easily move backwards
  if (detectorid == DetId::HGCalEE)
  {
    //First layer as no previous cells, point to 0.
    if (layerid == 1)
    {
      LOG(DEBUG) << "A"
                 << "\n";
      return DetId(0);
    }
    else
    {
      LOG(DEBUG) << "B"
                 << "\n";
      // Get the cell with the least x,y distance in the n-1th layer of the EE.
      return findCellCloseToXYpos(cellId, DetId::HGCalEE, subdetid, layerid - 1).first;
    }
  }
  if (detectorid == DetId::HGCalHSi || detectorid == DetId::HGCalHSc)
  {
    // from the frist layer in the hadronic part go back to the last layer in the EE
    if (layerid == 1)
    {
      LOG(DEBUG) << "D"
                 << "\n";
      return findCellCloseToXYpos(cellId, DetId::HGCalEE, subdetid, 28).first;
    }
    // HSc starts with layer 9, so for all cells from layer <10 we can just search in the Si part
    else if (layerid < 10)
    {
      LOG(DEBUG) << "C"
                 << "\n";
      return findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, layerid - 1).first;
    }

    // for layer > 10, both detectors need to be searched for the closest cells
    else
    {
      LOG(DEBUG) << "E"
                 << "\n";
      auto [sicanid, deltasi] = findCellCloseToXYpos(cellId, DetId::HGCalHSi, subdetid, layerid - 1);
      auto [sccanid, deltasc] = findCellCloseToXYpos(cellId, DetId::HGCalHSc, subdetid, layerid - 1);
      if (deltasi < deltasc)
      {
        return sicanid;
      }
      else
      {
        return sccanid;
      }
    }
  }
  // LOG(ERROR)  << hash << ""<<"\n"; //doenst work ?
  std::cout << hash;
  LOG(ERROR) << "Wont find neighbor. This part should never be reached."
             << "\n";
  exit(EXIT_FAILURE);
  return DetId(0);
}

// This is the function that does the search within the given detector/layer
std::pair<DetId, float> GeoExtractor::findCellCloseToXYpos(
    DetId originCellDetID,
    unsigned int targetdetectorid,
    unsigned int targetsubdetid,
    unsigned int targetlayerid)
{
  LOG(DEBUG) << "Start findCellCloseToXYpos for" << originCellDetID.rawId() << "\n";

  //Get xy from the origin cell
  float x = recHitTools.getPosition(originCellDetID).x();
  float y = recHitTools.getPosition(originCellDetID).y();
  std::pair<int, int> targetwaferortileid;
  std::pair<int, int> targetcellid;

  // get the position of the cell, that we are searching the neighbor for.
  if (isSiliconDet(originCellDetID.det()))
  {
    LOG(DEBUG) << "findCellCloseToXYpos:SiliconPos"
               << "\n";
    targetwaferortileid = recHitTools.getWafer(originCellDetID);
    targetcellid = recHitTools.getCell(originCellDetID);
  }
  else
  {
    LOG(DEBUG) << "findCellCloseToXYpos:Scintillator"
               << "\n";
    HGCScintillatorDetId scid = HGCScintillatorDetId(originCellDetID);
    targetwaferortileid = scid.ietaphi();
    targetcellid = std::make_pair(0, 0);
  }

  LOG(DEBUG) << "search from cell"
             << "\n";
  LOG(DEBUG) << *getCellptr(originCellDetID) << "\n";

  // The cell that is the closest, to be replace by close cells
  // if getcell can't find the wafer / cell it starts with the (0,0) coordinates
  DetId closest_cellDetId = getstartcell(targetdetectorid, targetsubdetid, targetlayerid, targetwaferortileid, targetcellid);
  Cell *closest_cellptr = getCellptr(closest_cellDetId);

  LOG(DEBUG) << "3: closest_cell"
             << "\n";

  float x_cur = closest_cellptr->x;
  float y_cur = closest_cellptr->y;
  float d_cur = (x_cur - x) * (x_cur - x) + (y_cur - y) * (y_cur - y);
  bool improvement = true;
  LOG(DEBUG) << "3.5"
             << "\n";
  while (improvement)
  {
    improvement = false;

    for (DetId &neighbor : closest_cellptr->neighbors)
    {
      LOG(DEBUG) << "    4"
                 << "\n";
      Cell *nptr = getCellptr(neighbor);
      LOG(DEBUG) << "\t" << getCellHashKeys(neighbor) << "\n";
      LOG(DEBUG) << "\t" << *nptr << "\n";
      float xn = nptr->x;
      float yn = nptr->y;
      float dn = (xn - x) * (xn - x) + (yn - y) * (yn - y);
      if (dn < d_cur)
      {
        LOG(DEBUG) << "\t"
                   << "improvement: " << dn << "<" << d_cur;
        improvement = true;
        closest_cellptr = nptr;
        d_cur = dn;
        break;
      }
    }
  }
  LOG(DEBUG) << "no improvement, exiting"
             << "\n";
  return std::make_pair(closest_cellptr->globalid, d_cur);
}

DetId GeoExtractor::getstartcell(unsigned int detectorid, unsigned int subdetid, unsigned int layerid, std::pair<int, int> waferortileid, std::pair<int, int> cellid)
{
  if (detcol.detectors.find(detectorid) == detcol.detectors.end())
  {
    LOG(ERROR) << "No such detector:\n"
               << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
    exit(EXIT_FAILURE);
  }
  Det &detector = detcol.detectors[detectorid];

  if (detector.subdetectors.find(subdetid) == detector.subdetectors.end())
  {
    LOG(ERROR) << "No such subdetector:\n"
               << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
    exit(EXIT_FAILURE);
  }
  Subdet &subdet = detector.subdetectors[subdetid];

  if (subdet.layers.find(layerid) == subdet.layers.end())
  {
    LOG(ERROR) << "No such layer:\n"
               << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
    exit(EXIT_FAILURE);
  }
  Layer &layer = subdet.layers[layerid];

  // for the EE and HSi part the cells are accessed via the wafers
  if (isSiliconDet(detectorid))
  {
    if (layer.wafers.begin() == layer.wafers.end())
    {
      LOG(ERROR) << "Wafer is empty:\n"
                 << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
      exit(EXIT_FAILURE);
    }

    //find the wafer, it may not exist in this layer, in this case search from the first in the list.
    if (layer.wafers.find(waferortileid) == layer.wafers.end())
    {
      LOG(DEBUG) << "Cant find wafer, using first wafer."
                 << "\n";
      waferortileid = layer.wafers.begin()->first;
    }
    Wafer &wafer = layer.wafers[waferortileid];

    //Find the cell
    if (wafer.cells.find(cellid) == wafer.cells.end())
    {
      LOG(DEBUG) << "Cant find cell, using first."
                 << "\n";
      cellid = wafer.cells.begin()->first;
    }
    Cell &cell = wafer.cells[cellid];

    LOG(DEBUG) << "search: " << detectorid << " " << subdetid << " " << layerid << " " << waferortileid << " " << cellid << "\n";

    LOG(DEBUG) << cell << "\n";
    HGCSiliconDetId siid = HGCSiliconDetId(cell.globalid);
    LOG(DEBUG) << "startcell: "
               << " " << cell.globalid.det() << " " << cell.globalid.subdetId() << " " << recHitTools.getLayer(cell.globalid) << " " << siid.waferUV() << " " << siid.cellUV() << "\n";

    return cell.globalid;
  }
  else
  {
    if (layer.tiles.begin() == layer.tiles.end())
    {
      LOG(ERROR) << "Tiles is empty\n"
                 << printcell(detectorid, subdetid, layerid, waferortileid, cellid);
      exit(EXIT_FAILURE);
    }
    if (layer.tiles.find(waferortileid) == layer.tiles.end())
    {
      waferortileid = layer.tiles.begin()->first;
    }
    return layer.tiles[waferortileid].globalid;
  }
  LOG(ERROR) << "Dead End1"
             << "\n";
  exit(EXIT_FAILURE);
  return DetId(0);
}
