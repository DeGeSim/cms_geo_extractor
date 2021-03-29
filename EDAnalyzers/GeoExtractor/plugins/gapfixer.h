#pragma once

// wraps findNextCell and loops over the ids
void GeoExtractor::fixGap(std::vector<DetId> &v_validHGCalIds)
{
  LOG(INFO) << "stating xposlist setup\n";
  setupXLists();
  LOG(INFO) << "xposlist setup done\n";
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

    altassingGapNeighbors(cellptr);

    DetId res;
    // //Continue to look for neighbors over the gap, until the number is sufficient or the Id0 is reached.
    // while (
    //     (iterId.det() == DetId::HGCalHSi && (int)cellptr->getAllNeighbors().size() < 6) ||
    //     (iterId.det() == DetId::HGCalHSc && (int)cellptr->getAllNeighbors().size() < 4))
    // {
    //   res = assingGapNeighbors(cellptr);
    //   if (res == DetId(0))
    //   {
    //     break;
    //   }

    //   cellptr->gapneighbors.insert(res);
    // }
  }
  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    //keep track of the number of added neighbors
    DetId iterId = v_validHGCalIds[i];
    Cell *cellptr = getCellPtr(iterId);
    treeOutput->ngapneighbors.push_back((unsigned int)cellptr->gapneighbors.size());
  }
}

DetId GeoExtractor::assingGapNeighbors(Cell *cellptr)
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
    LOG(ERROR) << "ERROR assingGapNeighbors: wrong detector type: " << detectorid;
    exit(EXIT_FAILURE);
  }

  auto [candidate, delta] = searchInLayer(cellptr->globalid, hash, targetdetectorid, subdetid, layerid, true);

  LOG(INFO) << "For origin " << cellptr->globalid.rawId() << ": ";
  LOG(INFO) << hash << "\n";
  LOG(INFO) << "For target " << candidate.rawId() << ": ";
  LOG(INFO) << getCellHash(candidate) << "\n";

  if (cellptr->gapneighbors.find(candidate) != cellptr->gapneighbors.end())
  {
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

void GeoExtractor::setupXLists()
{
  for (auto &[detectorid, det] : detcol.detectors)
  {
    xdistmap[detectorid];
    if (detectorid == DetId::HGCalEE)
    {
      continue;
    }
    for (auto &[subdetid, subdet] : det.subdetectors)
    {
      xdistmap[detectorid][subdetid];
      for (auto &[layerid, layer] : subdet.layers)
      {
        if (layerid < 9)
        {
          continue;
        }
        xdistmap[detectorid][subdetid][layerid];
        if (detectorid == DetId::HGCalHSi)
        {
          LOG(DEBUG) << detectorid << layerid << "wafers:" << (int)layer.wafers.size() << "\n";
          for (auto &[waferid, wafer] : layer.wafers)
          {
            LOG(DEBUG) << detectorid << layerid << waferid << "cells:" << (int)wafer.cells.size() << "\n";
            for (auto &[cellid, cell] : wafer.cells)
            {
              xdistmap[detectorid][subdetid][layerid].push_back(std::make_tuple(cell.x, &cell));
            }
          }
        }
        else
        {
          LOG(DEBUG) << detectorid << layerid << "tiles:" << (int)layer.tiles.size() << "\n";
          for (auto &[cellid, cell] : layer.tiles)
          {
            xdistmap[detectorid][subdetid][layerid].push_back(std::make_tuple(cell.x, &cell));
          }
        }
        std::vector<PosListTup> &xL = xdistmap[detectorid][subdetid][layerid];
        std::sort(xL.begin(), xL.end());
      }
    }
  }
}

void GeoExtractor::altassingGapNeighbors(Cell *cellptr)
{
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = getCellHash(cellptr->globalid);

  //The ist of the x values is presorted, the higher the index, the higher the number.
  std::vector<PosListTup> &xL = xdistmap[detectorid][subdetid][layerid];

  //These are the boundaries with cells, that could fulfill the distance condition,
  //they are initialized with the first and the last element
  std::vector<PosListTup>::iterator upper = xL.end() - 1;
  std::vector<PosListTup>::iterator lower = xL.begin();

  // Not find any element, that fullfils that condition in a binary search.
  bool cond = true;
  std::vector<PosListTup>::iterator current;
  double xdelta;
  bool inUpperRange, inLowerRange, higher;
  while (cond)
  {
    current = upper + (int)(std::distance(upper, lower) / 2);
    LOG(DEBUG) << "Current pos" << *current << "\n";
    auto &[curx, curptr] = *current;
    xdelta = std::abs(curx - cellptr->x);
    higher = curptr->x > cellptr->x;
    inUpperRange = higher && (xdelta < maxDeltaHScHSiGap);
    inLowerRange = !higher && (xdelta < maxDeltaHScHSiGap);

    //Stop if we have found an element in the range
    if (xdelta < maxDeltaHScHSiGap)
    {
      cond = false;
      LOG(DEBUG) << "Done! targetx:" << cellptr->x << " current:" << *current;
      LOG(DEBUG) << " upper:" << *upper;
      LOG(DEBUG) << " lower:" << *lower << "\n";
    }
    //if current is not with the lower boundery, set the lower boundery to current
    else if (inUpperRange)
    {
      LOG(DEBUG) << " A ";
      lower = current;
    }
    //vice versa
    else if (inLowerRange)
    {
      LOG(DEBUG) << " B ";
      upper = current;
    }
    else
    {
      //If current is neither in the lower nor in the higher boundery, but current>cell
      //then set the lower boundery to current
      if (higher)
      {
        LOG(DEBUG) << " C ";
        upper = current;
      }
      //vice versa
      else
      {
        LOG(DEBUG) << " D ";
        lower = current;
      }
    }
  }

  // Now, lower points to an element under the range (or on the edge)
  // upper points to an element over the range (or on the edge)
  // current is within the range
  LOG(DEBUG) << "length range " << std::distance(upper, lower) << "\n";
  // std::vector<PosListTup>::iterator lowercp = lower;
  // std::vector<PosListTup>::iterator uppercp = upper;
  // int raiselower = 0;
  // int lowerupper = 0;

  // while (!rangecond(lower, cellptr))
  // {
  //   raiselower++;
  //   lower++;
  // }
  // while (!rangecond(upper, cellptr))
  // {
  //   lowerupper++;
  //   upper--;
  // }
  // LOG(DEBUG) << "raiselower " << raiselower << " lowerupper:" << lowerupper << "\n";

  // LOG(DEBUG) << "Compare lowerbound: " << std::distance(lower, lower) << "\n";
  // LOG(DEBUG) << "Compare upperbound: " << std::distance(upper, upper) << "\n";

  std::vector<PosListTup>::iterator inner = current;
  std::vector<PosListTup>::iterator middle;
  int iterations = 0;
  //Binary search to put lower to the lower edge of the boundery.
  while (!rangecond(lower, cellptr))
  {
    middle = lower + std::distance(lower, inner) / 2;
    LOG(DEBUG) << "l " << (lower - xL.begin()) << " m " << (middle - xL.begin()) << " i " << (inner - xL.begin()) << "\n";
    if (!rangecond(middle, cellptr))
    {
      lower = middle;
    }
    else
    {
      inner = middle;
    }
    if (lower == middle && lower + 1 == inner)
    {
      lower = inner;
      break;
    }
    iterations++;
    if (iterations > 14)
    {
      exit(0);
    }
  }
  LOG(DEBUG) << "Lower bound found in " << iterations << " iterations.\n";

  inner = current;
  iterations = 0;
  //Binary search to put upper to the upper edge of the boundery.
  while (!rangecond(upper, cellptr))
  {
    //The difference is negative
    middle = upper + std::distance(upper, inner) / 2;
    LOG(DEBUG) << "u " << (upper - xL.begin()) << " m " << (middle - xL.begin()) << " i " << (inner - xL.begin()) << "\n";
    if (!rangecond(middle, cellptr))
    {
      upper = middle;
    }
    else
    {
      inner = middle;
    }
    if (upper == middle && upper - 1 == inner)
    {
      upper = inner;
      LOG(DEBUG) << "u " << (upper - xL.begin()) << " m " << (middle - xL.begin()) << " i " << (inner - xL.begin()) << "\n";
      break;
    }
    iterations++;
    if (iterations > 14)
    {
      exit(0);
    }
  }
  LOG(DEBUG) << "Upper bound found in " << iterations << "iterations.\n";

  LOG(DEBUG) << "distance lower upper " << std::distance(lower, upper) << "\n";
  if (std::distance(lower, lower) != 0 || std::distance(upper, upper))
  {
    exit(1);
  }
  //Now upper points to the highest element in the range, lower to the lowerest

  //Now also check the condition in y and save the cellpointer to a vector
  std::vector<PosListTup> v_newGapNeighbors;
  for (vector<PosListTup>::iterator candit = lower; candit <= upper; candit++)
  {
    Cell *cellcanditptr = std::get<1>(*candit);
    double dx = (cellcanditptr->x - cellptr->x);
    double dy = (cellcanditptr->y - cellptr->y);
    double delta = dx * dx + dy * dy;
    if (delta < maxDeltaHScHSiGap)
    {
      v_newGapNeighbors.push_back(std::make_tuple(delta, cellcanditptr));
    }
  }

  LOG(INFO) << "number of new neighbors:" << (int)v_newGapNeighbors.size() << "\n";
  std::sort(v_newGapNeighbors.begin(), v_newGapNeighbors.end());

  int maxneighbors = (isSiliconDet(cellptr->globalid.det())) ? 8 : 6;
  int curneighbors;
  int iadded = 0;
  for (auto &[delta, gapneighborptr] : v_newGapNeighbors)
  {
    curneighbors = (int)cellptr->neighbors.size() + (int)cellptr->gapneighbors.size();
    if (curneighbors >= maxneighbors)
    {
      LOG(INFO) << "Stop adding neigbors, because total number is " << curneighbors << "\n";
      break;
    }
    cellptr->gapneighbors.insert(gapneighborptr->globalid);
    iadded++;
    LOG(INFO) << cellptr->globalid.rawId() << ": Adding (" << iadded << ") " << gapneighborptr->globalid.rawId() << "\n";
  }
}

bool GeoExtractor::rangecond(std::vector<PosListTup>::iterator iter, Cell *cellptr)
{
  return (std::abs(std::get<0>(*iter) - cellptr->x) < maxDeltaHScHSiGap);
}