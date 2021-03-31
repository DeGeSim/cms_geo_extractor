#pragma once

// wraps findNextCell and loops over the ids
void GeoExtractor::fixGap(std::vector<DetId> &v_validHGCalIds)
{
  LOG(INFO) << "stating xposlist setup\n";
  setupXLists();
  LOG(INFO) << "xposlist setup done\n";
  LOGCFG.level = DEBUG;
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
    if (iterId.det() == DetId::HGCalHSi)
    {
      if ((int)cellptr->getAllNeighbors().size() >= simaxneighbors)
        continue;
    }
    if (iterId.det() == DetId::HGCalHSc)
    {
      if ((int)cellptr->getAllNeighbors().size() >= scmaxneighbors)
        continue;
    }

    LOG(DEBUG) << "Assigning gapneighbors for " << cellptr->globalid.rawId() << " (#" << cellptr->getAllNeighbors().size() << "): \n";
    LOG(DEBUG) << *cellptr << "\n";

    assingGapNeighbors(cellptr);

    DetId res;
  }
  LOGCFG.level = internalDebugLevel;
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
        std::sort(
            xdistmap[detectorid][subdetid][layerid].begin(),
            xdistmap[detectorid][subdetid][layerid].end());
      }
    }
  }
}

void GeoExtractor::assingGapNeighbors(Cell *cellptr)
{
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = getCellHash(cellptr->globalid);

  LOG(DEBUG) << "assingGapNeighbors cell " << cellptr->globalid.rawId() << "\n";
  LOG(DEBUG) << getCellHash(cellptr->globalid) << "\n";
  if (detectorid == DetId::HGCalHSi)
  {
    //skip if layer doesnt exist or has not been initialized because it is empty.
    std::map<int, Layer> &layers = detcol.detectors[DetId::HGCalHSc].subdetectors[subdetid].layers;
    if (layers.find(layerid) == layers.end())
      return;
  }
  if (detectorid == DetId::HGCalHSc)
  {
    //skip if layer doesnt exist or has not been initialized because it is empty.
    std::map<int, Layer> &layers = detcol.detectors[DetId::HGCalHSc].subdetectors[subdetid].layers;
    if (layers.find(layerid) == layers.end())
      return;
  }

  //xL is the of the (x,cellptr) pair presorted by x, the higher the index, the higher the number.
  std::vector<PosListTup> xL = (detectorid == DetId::HGCalHSi) ? xdistmap[DetId::HGCalHSc][subdetid][layerid] : xdistmap[DetId::HGCalHSi][subdetid][layerid];

  //These are the boundaries with cells, that could fulfill the distance condition,
  //they are initialized with the first and the last element
  std::vector<PosListTup>::iterator upper = xL.end() - 1;
  std::vector<PosListTup>::iterator lower = xL.begin();

  if (std::get<0>(*upper) < cellptr->x - maxDeltaHScHSiGap)
  {
    LOG(DEBUG) << "xL highest element xpos is less than the target xpos - maxDeltaHScHSiGap, aborting gapneighbor search for this cell.\n";
    return;
  }
  if (std::get<0>(*lower) > cellptr->x + maxDeltaHScHSiGap)
  {
    LOG(DEBUG) << "xL lowest element xpos is lower than the target xpos + maxDeltaHScHSiGap, aborting gapneighbor search for this cell.\n";
    return;
  }

  ///////////////////
  // Now find any element in the range xdelta around the xposition, a binary search.
  ///////////////////

  LOG(DEBUG) << "Searching between " << *lower << " and " << *upper << "\n";

  bool cond = true;
  std::vector<PosListTup>::iterator current;
  double xdelta;
  bool inUpperRange, inLowerRange, higher;
  int ifoo = 0;
  while (cond)
  {
    ifoo++;
    if (ifoo > 30)
      exit(EXIT_FAILURE);
    current = upper + (int)(std::distance(upper, lower) / 2);

    // LOG(DEBUG) << "Current pos" << *current << "\n";
    LOG(DEBUG) << "l " << (lower - xL.begin()) << " c " << (current - xL.begin()) << " u " << (upper - xL.begin()) << "\n";

    auto &[curx, curptr] = *current;
    xdelta = std::abs(curx - cellptr->x);
    higher = curptr->x > cellptr->x;
    inUpperRange = higher && (xdelta <= maxDeltaHScHSiGap);
    inLowerRange = !higher && (xdelta <= maxDeltaHScHSiGap);

    //Stop if we have found an element in the range
    if (xdelta <= maxDeltaHScHSiGap)
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
        //Special Case handling if only the element on the lower end of the list is in range.
        if ((upper == current && upper - 1 == lower) && (lower - xL.begin() == 0))
        {
          LOG(DEBUG) << " C1 ";
          current = lower;
          break;
        }
        else
        {
          upper = current;
        }
      }
      //vice versa
      else
      {
        LOG(DEBUG) << " D ";
        //Special Case handling if only the element on the upper end of the list is in range.
        if ((lower == current && upper - 1 == lower) && (xL.end() - upper == 1))
        {
          LOG(DEBUG) << " D1 ";
          current = upper;
          break;
        }
        else
        {
          lower = current;
        }
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
  LOG(DEBUG) << "For " << cellptr->globalid.rawId() << "\n";
  for (vector<PosListTup>::iterator candit = lower; candit <= upper; candit++)
  {
    Cell *cellcanditptr = std::get<1>(*candit);
    if (cellcanditptr->globalid == cellptr->globalid)
      continue;

    double dx = (cellcanditptr->x - cellptr->x);
    double dy = (cellcanditptr->y - cellptr->y);
    double delta = std::sqrt(dx * dx + dy * dy);
    if (delta < maxDeltaHScHSiGap)
    {
      LOG(DEBUG) << delta;
      LOG(DEBUG) << " from (" << cellptr->x << "," << cellptr->y;
      LOG(DEBUG) << ") vs (" << cellcanditptr->x << "," << cellcanditptr->y << ")\n";
      v_newGapNeighbors.push_back(std::make_tuple(delta, cellcanditptr));
    }
  }
  LOG(DEBUG) << "\n";

  LOG(DEBUG) << "number of opetential new neighbors:" << (int)v_newGapNeighbors.size() << "\n";
  std::sort(v_newGapNeighbors.begin(), v_newGapNeighbors.end());

  int maxneighbors = (isSiliconDet(cellptr->globalid.det())) ? simaxneighbors : scmaxneighbors;
  int curneighbors;
  int iadded = 0;
  for (auto &[delta, gapneighborptr] : v_newGapNeighbors)
  {
    curneighbors = (int)cellptr->neighbors.size() + (int)cellptr->gapneighbors.size();
    if (curneighbors >= maxneighbors)
    {

      LOG(DEBUG) << "Stop adding neigbors, because total number is " << curneighbors << "\n";
      break;
    }

    LOG(DEBUG) << cellptr->globalid.rawId() << ": Adding (" << iadded << ") " << gapneighborptr->globalid.rawId() << " delta " << delta << "\n";

    LOG(DEBUG) << cellptr->globalid.rawId() << " " << (int)cellptr->getAllNeighbors().size() << " " << (int)cellptr->neighbors.size() << " " << (int)cellptr->gapneighbors.size() << " "
               << "\n";

    cellptr->gapneighbors.insert(gapneighborptr->globalid);
    iadded++;
  }
}

bool GeoExtractor::rangecond(std::vector<PosListTup>::iterator iter, Cell *cellptr)
{
  return (std::abs(std::get<0>(*iter) - cellptr->x) < maxDeltaHScHSiGap);
}