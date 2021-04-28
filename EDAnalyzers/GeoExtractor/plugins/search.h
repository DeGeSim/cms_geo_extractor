#pragma once

// wraps findNextCell and loops over the ids
void GeoExtractor::assignZNeighbors(std::vector<DetId> &v_validHGCalIds)
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
      LOG(DEBUG) << "Assigning z neighbors " << i << "\n";
    }
    Cell *cellptr = getCellPtr(iterId);
    validateId(cellptr->globalid);
    cellptr->next = findNextCell(cellptr->globalid);
    cellptr->previous = findPreviousCell(cellptr->globalid);
  }
}

// This function decices in which layer in which detector(s) to search for the neigbor.
// The the act lifting is done by the searchInLayer
DetId GeoExtractor::findNextCell(DetId cellId)
{
  LOG(DEBUG) << "Start findNextCell\n";
  CellHash hash = getCellHash(cellId);
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
  LOG(DEBUG) << "find cell for id" << cellId.rawId() << "\n";
  LOG(DEBUG) << getCellHash(cellId) << "\n";
  //
  // HGCalEE = 8, layer 1-28
  // HGCalHSi = 9, layer 1-22
  // HGCalHSc = 10, layer 9-22
  // HGCalTrigger = 11 X
  
  int direction = recHitTools.getPosition(cellId).z() > 0 ? 1 : -1;

  // For the ee cal we can easily move forward
  if (detectorid == DetId::HGCalEE)
  {
    if (std::abs(layerid) < 28)
    {
      LOG(DEBUG) << "A\n";
      return searchInLayer(cellId, hash, detectorid, subdetid, layerid + 1*direction).first;
    }
    // to the the next detector EE-> HSi
    else
    {
      LOG(DEBUG) << "B\n";
      return searchInLayer(cellId, hash, DetId::HGCalHSi, subdetid, 1*direction).first;
    }
  }
  //For the hadronic part
  if (detectorid == DetId::HGCalHSi || detectorid == DetId::HGCalHSc)
  {
    // for layer <8 all cells we can just search in the Si part
    if (std::abs(layerid) < 8)
    {
      LOG(DEBUG) << "C\n";
      return searchInLayer(cellId, hash, DetId::HGCalHSi, subdetid, layerid + 1*direction).first;
    }
    // return 0 in the last layer
    else if (std::abs(layerid) == 22)
    {
      LOG(DEBUG) << "D\n";
      return DetId(0);
    }
    else
    {
      LOG(DEBUG) << "E1\n";

      // Avoid searching bot subdetectors if one of them doesnt have cells in the layer:
      Det &detectorsi = detcol.detectors[DetId::HGCalHSi];
      Subdet &subdetsi = detectorsi.subdetectors[subdetid];
      if (subdetsi.layers.find(layerid + 1*direction) == subdetsi.layers.end())
      {
        LOG(DEBUG) << "E1-1\n";
        return searchInLayer(cellId, hash, DetId::HGCalHSc, subdetid, std::abs(layerid) + 1*direction).first;
      }
      LOG(DEBUG) << "E2\n";
      Det &detectorsc = detcol.detectors[DetId::HGCalHSc];
      Subdet &subdetsc = detectorsc.subdetectors[subdetid];
      if (subdetsc.layers.find(layerid + 1*direction) == subdetsc.layers.end())
      {
        LOG(DEBUG) << "E1-1\n";
        return searchInLayer(cellId, hash, DetId::HGCalHSi, subdetid, layerid + 1*direction).first;
      }

      LOG(DEBUG) << "E3\n";
      auto [sicanid, deltasi] = searchInLayer(cellId, hash, DetId::HGCalHSi, subdetid, layerid + 1*direction);
      LOG(DEBUG) << "E4\n";
      auto [sccanid, deltasc] = searchInLayer(cellId, hash, DetId::HGCalHSc, subdetid, layerid + 1*direction);

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
  LOG(ERROR) << "Wont find neighbor. This part should never be reached.\n";
  exit(EXIT_FAILURE);
  return DetId(0);
}

DetId GeoExtractor::findPreviousCell(DetId cellId)
{
  LOG(DEBUG) << "Start findPreviousCell\n";
  CellHash hash = getCellHash(cellId);
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
  LOG(DEBUG) << "find cell for id" << cellId.rawId() << "\n";
  LOG(DEBUG) << hash << "\n";
  //
  // HGCalEE = 8, layer 1-28
  // HGCalHSi = 9, layer 1-22
  // HGCalHSc = 10, layer 9-22
  // HGCalTrigger = 11 X
  int direction = (layerid>0) ? 1:-1;

  // For the ee cal we can easily move backwards
  if (detectorid == DetId::HGCalEE)
  {
    //First layer as no previous cells, point to 0.
    if (std::abs(layerid) == 1)
    {
      LOG(DEBUG) << "A\n";
      return DetId(0);
    }
    else
    {
      LOG(DEBUG) << "B\n";
      // Get the cell with the least x,y distance in the n-1th layer of the EE.
      return searchInLayer(cellId, hash, DetId::HGCalEE, subdetid, layerid - direction).first;
    }
  }
  if (detectorid == DetId::HGCalHSi || detectorid == DetId::HGCalHSc)
  {
    // from the first layer in the hadronic part go back to the last layer in the EE
    if (std::abs(layerid) == 1)
    {
      LOG(DEBUG) << "D\n";
      return searchInLayer(cellId, hash, DetId::HGCalEE, subdetid, direction*28).first;
    }
    // HSc starts with layer 9, so for all cells from layer <10 we can just search in the Si part
    else if (std::abs(layerid) < 10)
    {
      LOG(DEBUG) << "C\n";
      return searchInLayer(cellId, hash, DetId::HGCalHSi, subdetid, layerid - direction).first;
    }

    // for layer > 10, both detectors need to be searched for the closest cells
    else
    {
      LOG(DEBUG) << "E\n";

      // Avoid searching bot subdetectors if one of them doenst have cells in the layer:
      Det &detectorsi = detcol.detectors[DetId::HGCalHSi];
      Subdet &subdetsi = detectorsi.subdetectors[subdetid];
      if (subdetsi.layers.find(layerid - direction) == subdetsi.layers.end())
      {
        return searchInLayer(cellId, hash, DetId::HGCalHSc, subdetid, layerid - direction).first;
      }
      Det &detectorsc = detcol.detectors[DetId::HGCalHSc];
      Subdet &subdetsc = detectorsc.subdetectors[subdetid];
      if (subdetsc.layers.find(layerid - direction) == subdetsc.layers.end())
      {
        return searchInLayer(cellId, hash, DetId::HGCalHSi, subdetid, layerid - direction).first;
      }

      auto [sicanid, deltasi] = searchInLayer(cellId, hash, DetId::HGCalHSi, subdetid, layerid - direction);
      auto [sccanid, deltasc] = searchInLayer(cellId, hash, DetId::HGCalHSc, subdetid, layerid - direction);
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
  // LOG(ERROR)  << hash << "\n"; //doenst work ?
  std::cout << hash;
  LOG(ERROR) << "Wont find neighbor. This part should never be reached.\n";
  exit(EXIT_FAILURE);
  return DetId(0);
}

// This is the function that does the search within the given detector/layer
std::pair<DetId, float> GeoExtractor::searchInLayer(
    DetId originCellDetID,
    CellHash hash,
    unsigned int targetdetectorid,
    unsigned int targetsubdetid,
    int targetlayerid)
{
  LOG(DEBUG) << "Start searchInLayer for Id " << originCellDetID.rawId() << "\n";

  std::pair<int, int> targetwaferortileid = std::get<3>(hash);
  std::pair<int, int> targetcellid = std::get<4>(hash);

  Cell *originCellptr = getCellPtr(originCellDetID);
  LOG(DEBUG) << "search from cell\n";
  LOG(DEBUG) << *originCellptr << "\n";

  // Prepare for the getStartCell function:
  // Create a hash with the desired target coordinates as an entry point.
  CellHash targethash = std::make_tuple(targetdetectorid, targetsubdetid, targetlayerid, targetwaferortileid, targetcellid);

  DetId closest_cellDetId = getStartCell(targethash);
  // if getStartCell can't find the wafer / cell it starts with the (u,v)=(0,0) coordinates
  Cell *closest_cellptr = getCellPtr(closest_cellDetId);

  LOG(DEBUG) << "3: closest_cell\n";
  //Get xy from the origin cell
  float x = recHitTools.getPosition(originCellDetID).x();
  float y = recHitTools.getPosition(originCellDetID).y();

  float x_cur = closest_cellptr->x;
  float y_cur = closest_cellptr->y;
  float d_cur = (x_cur - x) * (x_cur - x) + (y_cur - y) * (y_cur - y);
  bool improvement = true;
  LOG(DEBUG) << "3.5\n";
  while (improvement)
  {
    improvement = false;

    for (auto neighbor : closest_cellptr->neighbors)
    {
      LOG(DEBUG) << "    4\n";
      Cell *nptr = getCellPtr(neighbor);
      LOG(DEBUG) << "\t" << getCellHash(neighbor) << "\n";
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
  LOG(DEBUG) << "no improvement, exiting\n";
  return std::make_pair(closest_cellptr->globalid, d_cur);
}

DetId GeoExtractor::getStartCell(CellHash hash)
{
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
  if (detcol.detectors.find(detectorid) == detcol.detectors.end())
  {
    LOG(ERROR) << "No such detector:\n";
    LOG(ERROR) << printCell(detectorid, subdetid, layerid, waferortileid, cellid);
    exit(EXIT_FAILURE);
  }
  Det &detector = detcol.detectors[detectorid];

  if (detector.subdetectors.find(subdetid) == detector.subdetectors.end())
  {
    LOG(ERROR) << "No such subdetector:\n";
    LOG(ERROR) << printCell(detectorid, subdetid, layerid, waferortileid, cellid);
    exit(EXIT_FAILURE);
  }
  Subdet &subdet = detector.subdetectors[subdetid];

  if (subdet.layers.find(layerid) == subdet.layers.end())
  {
    LOG(WARN) << "No such layer:\n";
    LOG(WARN) << printCell(detectorid, subdetid, layerid, waferortileid, cellid);
    exit(EXIT_FAILURE);
  }
  Layer &layer = subdet.layers[layerid];

  // for the EE and HSi part the cells are accessed via the wafers
  if (isSiliconDet(detectorid))
  {
    if (layer.wafers.begin() == layer.wafers.end())
    {
      LOG(ERROR) << "Wafer is empty:\n";
      LOG(ERROR) << printCell(detectorid, subdetid, layerid, waferortileid, cellid);
      exit(EXIT_FAILURE);
    }

    //find the wafer, it may not exist in this layer, in this case search from the first in the list.
    if (layer.wafers.find(waferortileid) == layer.wafers.end())
    {
      LOG(DEBUG) << "Cant find wafer, using first wafer.\n";
      waferortileid = layer.wafers.begin()->first;
    }
    Wafer &wafer = layer.wafers[waferortileid];

    //Find the cell
    if (wafer.cells.find(cellid) == wafer.cells.end())
    {
      LOG(DEBUG) << "Cant find cell, using first.\n";
      cellid = wafer.cells.begin()->first;
    }
    Cell *cellptr = &wafer.cells[cellid];

    LOG(DEBUG) << "search: \n";
    LOG(DEBUG) << printCell(detectorid, subdetid, layerid, waferortileid, cellid);

    LOG(DEBUG) << *cellptr << "\n";
    HGCSiliconDetId siid = HGCSiliconDetId(cellptr->globalid);
    LOG(DEBUG) << "startcell: \n";
    LOG(DEBUG) << printCell(cellptr->globalid.det(), cellptr->globalid.subdetId(), realLayerFromId(cellptr->globalid), siid.waferUV(), siid.cellUV());

    return cellptr->globalid;
  }
  else
  {
    if (layer.tiles.begin() == layer.tiles.end())
    {
      LOG(ERROR) << "Tiles is empty\n";
      LOG(ERROR) << printCell(detectorid, subdetid, layerid, waferortileid, cellid);
      exit(EXIT_FAILURE);
    }
    if (layer.tiles.find(waferortileid) == layer.tiles.end())
    {
      waferortileid = layer.tiles.begin()->first;
    }
    return layer.tiles[waferortileid].globalid;
  }
  LOG(ERROR) << "Dead End1\n";
  exit(EXIT_FAILURE);
  return DetId(0);
}
