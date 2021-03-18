#pragma once

// wraps findNextCell and loops over the ids
void GeoExtractor::fixNeighborsBoundery(std::vector<DetId> &v_validHGCalIds)
{
  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    //skip for HGCalEE and for layers where HSc and HSi dont overlap
    DetId iterId = v_validHGCalIds[i];

    if (iterId.det() == DetId::HGCalEE || recHitTools.getLayer(iterId) < 9)
    {
      continue;
    }
    Cell *cellptr = getCellPtr(iterId);
    //Skip if there is allready a sufficient number of neighbors
    if (iterId.det() == DetId::HGCalHSi && (int)cellptr->neighbors.size() >= 6)
    {
      continue;
    }
    if (iterId.det() == DetId::HGCalHSc && (int)cellptr->neighbors.size() >= 4)
    {
      continue;
    }

    DetId res;
    while ((iterId.det() == DetId::HGCalHSi && (int)cellptr->neighbors.size() < 6) || (iterId.det() == DetId::HGCalHSc && (int)cellptr->neighbors.size() < 4))
    {
      res = findGapNeighbors(cellptr);
      if (res == DetId(0))
      {
        break;
      }
      LOG(INFO) << "Adding neighbor for " << cellptr->globalid.rawId() << " (#" << cellptr->neighbors.size() << ") : " << res.rawId() << "\n";
      cellptr->neighbors.push_back(res);
    }
  }
}
DetId GeoExtractor::findGapNeighbors(Cell *cellptr)
{
  CellHash hash = getCellHash(cellptr->globalid);
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = hash;
  int targetdetectorid;
  if (detectorid == DetId::HGCalHSi)
  {
    targetdetectorid = DetId::HGCalHSc;
  }
  else
  {
    targetdetectorid = DetId::HGCalHSi;
  }
  auto [candidate, delta] = searchInLayer(cellptr->globalid, hash, targetdetectorid, subdetid, layerid, true);
  if (delta > maxDeltaHScHSiGap)
  {
    return DetId(0);
  }
  else
  {
    return candidate;
  }
}
