#pragma once

// wraps findNextCell and loops over the ids
void GeoExtractor::fixNeighborsBoundery(std::vector<DetId> &v_validHGCalIds)
{
  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    //skip for HGCalEE and for layers where HSc and HSi dont overlap
    DetId iterId = v_validHGCalIds[i];

    if (iterId.det() != DetId::HGCalHSi && iterId.det() != DetId::HGCalHSc)
    {
      continue;
    }
    if (recHitTools.getLayer(iterId) < 9)
    {
      continue;
    }
    Cell *cellptr = getCellPtr(iterId);
    //Skip if there is allready a sufficient number of neighbors
    if (iterId.det() == DetId::HGCalHSi && (int)cellptr->getAllNeighbors().size() >= 6)
    {
      continue;
    }
    if (iterId.det() == DetId::HGCalHSc && (int)cellptr->getAllNeighbors().size() >= 4)
    {
      continue;
    }
    // LOG(INFO) << "Detector " << iterId.det() << "\n";
    // LOG(INFO) << "Other Detector " << std::get<0>(getCellHash(cellptr->globalid)) << "\n";
    LOG(INFO) << "Adding neighbor for " << cellptr->globalid.rawId() << " (#" << cellptr->getAllNeighbors().size() << "): \n";
    LOG(INFO) << *cellptr << "\n";

    DetId res;
    //Continue to look for neighbors over the gap, until the number is sufficient or the Id0 is reached.
    while (
        (iterId.det() == DetId::HGCalHSi && (int)cellptr->getAllNeighbors().size() < 6) ||
        (iterId.det() == DetId::HGCalHSc && (int)cellptr->getAllNeighbors().size() < 4))
    {
      res = findGapNeighbors(cellptr);
      if (res == DetId(0))
      {
        break;
      }

      cellptr->gapneighbors.insert(res);
    }
  }
  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    //keep track of the number of added neighbors
    DetId iterId = v_validHGCalIds[i];
    Cell *cellptr = getCellPtr(iterId);
    treeOutput->ngapneighbors.push_back((unsigned int)cellptr->gapneighbors.size());
  }
}
DetId GeoExtractor::findGapNeighbors(Cell *cellptr)
{
  CellHash hash = getCellHash(cellptr->globalid);
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
  int targetdetectorid;
  switch (detectorid)
  {
  case DetId::HGCalHSi:
    targetdetectorid = DetId::HGCalHSc;
    break;
  case DetId::HGCalHSc:
    targetdetectorid = DetId::HGCalHSi;
    break;
  default:
    LOG(ERROR) << "ERROR findGapNeighbors: wrong detector type: " << detectorid;
    exit(EXIT_FAILURE);
  }

  LOGCFG.level = DEBUG;
  auto [candidate, delta] = searchInLayer(cellptr->globalid, hash, targetdetectorid, subdetid, layerid, true);
  LOGCFG.level = internalDebugLevel;
  
  LOG(INFO) << "For origin " << cellptr->globalid.rawId() <<": ";
  LOG(INFO) << hash << "\n";
  LOG(INFO) << "For target " << candidate.rawId() << ": ";
  LOG(INFO) << getCellHash(candidate) << "\n";
  
  if (cellptr->gapneighbors.find(candidate)!=cellptr->gapneighbors.end()) {
    LOG(INFO) << "Candidate is already part of the neighbors for this cell, stoping search.\n";
    return DetId(0);
  }
  else if (delta > maxDeltaHScHSiGap)
  {
    LOG(INFO) << "Gap " << delta << " vs " << maxDeltaHScHSiGap;
    LOG(INFO) << " => skipping.\n\n";
    return DetId(0);
  }
  else
  {
    LOG(INFO) << "Gap " << delta << " vs " << maxDeltaHScHSiGap;
    LOG(INFO) << " => adding.\n\n";
    return candidate;
  }
}
