#pragma once

// wraps findNextCell and loops over the ids
void GeoExtractor::fixGap(std::vector<DetId> &v_validHGCalIds)
{
  // LOGCFG.level = DEBUG;
  LOG(INFO) << "starting xposlist setup\n";
  setupXLists();
  LOG(INFO) << "xposlist setup done\n";
  for (int i = 0; i < (int)v_validHGCalIds.size(); i++)
  {
    //skip for HGCalEE and for layers where HSc and HSi dont overlap
    DetId iterId = v_validHGCalIds[i];

    //Cell with super irregular geometry
    // if (iterId.rawId() != 2752652289) continue;

    if (iterId.det() != DetId::HGCalHSc)
      continue;
    //Use the absolute of the layer for this case
    if (recHitTools.getLayer(iterId) < 9)
      continue;

    Cell *cellptr = getCellPtr(iterId);
    
    LOG(DEBUG) << "Assigning gapneighbors for " << cellptr->globalid.rawId() << " (#" << cellptr->getAllNeighbors().size() << "): \n";
    LOG(DEBUG) << *cellptr << "\n";

    assignGapNeighbors(cellptr);
  }
  // LOGCFG.level = internalDebugLevel;
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
        if (std::abs(layerid) < 9)
        {
          continue;
        }
        xdistmap[detectorid][subdetid][layerid];
        if (detectorid == DetId::HGCalHSi)
        {
          LOG(DEBUG) << "detectorid  " << detectorid << " layerid " << layerid << "wafers:" << (int)layer.wafers.size() << "\n";
          for (auto &[waferid, wafer] : layer.wafers)
          {
            LOG(DEBUG) << "detectorid  " << detectorid << " layerid " << layerid << "waferid " << waferid << "cells: " << (int)wafer.cells.size() << "\n";
            for (auto &[cellid, cell] : wafer.cells)
            {
              xdistmap[detectorid][subdetid][layerid].push_back(std::make_tuple(cell.x, &cell));
            }
          }
        }
        else
        {
          LOG(DEBUG) << "detectorid  " << detectorid << " layerid " << layerid << "tiles:" << (int)layer.tiles.size() << "\n";
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

void GeoExtractor::assignGapNeighbors(Cell *cellptr)
{
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = getCellHash(cellptr->globalid);

  LOG(DEBUG) << "GeoExtractor::assignGapNeighbors cell " << cellptr->globalid.rawId() << "\n";
  LOG(DEBUG) << *cellptr << "\n";
  LOG(DEBUG) << getCellHash(cellptr->globalid) << "\n";
  // Only cells from HGCalHSc are used now
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
  LOG(DEBUG) << "Skipcheck done.\n";
  //xL is the of the (x,cellptr) pair presorted by x, the higher the index, the higher x.
  std::vector<PosListTup> xL = (detectorid == DetId::HGCalHSi) ? xdistmap[DetId::HGCalHSc][subdetid][layerid] : xdistmap[DetId::HGCalHSi][subdetid][layerid];
  LOG(DEBUG) << "Got  xdistmap with size " << (int)xL.size() << "\n";

  // If the window is very small, the list may be empty.
  if ((int)xL.size()==0) return;
  //These are the boundaries with cells, that could fulfill the distance condition,
  //they are initialized with the first and the last element
  std::vector<PosListTup>::iterator upper = xL.end() - 1;
  std::vector<PosListTup>::iterator lower = xL.begin();

  //Abort if the x position of the cell can not be reached from the highest/lowest element in the list
  if (std::get<0>(*upper) < cellptr->x - maxSearchDelta)
  {
    LOG(DEBUG) << "xL highest element xpos is less than the target xpos - maxSearchDelta, aborting gapneighbor search for this cell.\n";
    return;
  }
  if (std::get<0>(*lower) > cellptr->x + maxSearchDelta)
  {
    LOG(DEBUG) << "xL lowest element xpos is lower than the target xpos + maxSearchDelta, aborting gapneighbor search for this cell.\n";
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
    inUpperRange = higher && (xdelta <= maxSearchDelta);
    inLowerRange = !higher && (xdelta <= maxSearchDelta);

    //Stop if we have found an element in the range
    if (xdelta <= maxSearchDelta)
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
  }
  LOG(DEBUG) << "Upper bound found in " << iterations << "iterations.\n";

  LOG(DEBUG) << "distance lower upper " << std::distance(lower, upper) << "\n";
  LOG(DEBUG) << "original cell: " << *cellptr << "\n";
  LOG(DEBUG) << "upper: " << *std::get<1>(*upper) << "\n";
  LOG(DEBUG) << "xposdiff: " << xposdiff(upper, cellptr) << "\n";
  LOG(DEBUG) << "lower: " << *std::get<1>(*lower) << "\n";
  LOG(DEBUG) << "xposdiff: " << xposdiff(lower, cellptr) << "\n";

  //Now upper points to the highest element in the range, lower to the lowerest
  //Now also check the condition in y and save the cellpointer to a vector
  std::vector<PosListTup> v_newGapNeighbors;
  LOG(DEBUG) << "After filtering for x, check all candidates for " << cellptr->globalid.rawId() << "\n";
  for (vector<PosListTup>::iterator candit = lower; candit <= upper; candit++)
  {
    Cell *cellcanditptr = std::get<1>(*candit);
    if (cellcanditptr->globalid == cellptr->globalid){
      LOG(ERROR) << "Error! cellcanditptr->globalid == cellptr->globalid\n";
      LOG(ERROR) << "cellcandit: \n" << *cellcanditptr;
      exit(EXIT_FAILURE);
    }

    //Filter for x and y now
    double dx = (cellcanditptr->x - cellptr->x);
    double dy = (cellcanditptr->y - cellptr->y);
    if (std::sqrt(dx * dx + dy * dy) > maxSearchDelta)
      continue;

    // Now check the real criterion (maxDeltaHScHSiGap)
    // with the minimum distance between the edges intead of the centers:
    double delta = cellsDelta(cellptr, cellcanditptr);
    if (delta < maxDeltaHScHSiGap)
    {
      LOG(DEBUG) << delta;
      LOG(DEBUG) << " from (" << cellptr->x << "," << cellptr->y;
      LOG(DEBUG) << ") vs (" << cellcanditptr->x << "," << cellcanditptr->y << ")\n";
      v_newGapNeighbors.push_back(std::make_tuple(delta, cellcanditptr));
    }
  }
  LOG(DEBUG) << "\n";

  LOG(DEBUG) << "number of potential new neighbors:" << (int)v_newGapNeighbors.size() << "\n";
  std::sort(v_newGapNeighbors.begin(), v_newGapNeighbors.end());

  int iadded = 0;
  for (auto &[delta, gapneighborptr] : v_newGapNeighbors)
  {
    LOG(DEBUG) << cellptr->globalid.rawId() << ": Adding (" << iadded << ") " << gapneighborptr->globalid.rawId() << " delta " << delta << "\n";

    LOG(DEBUG) << cellptr->globalid.rawId() << " " << (int)cellptr->getAllNeighbors().size() << " " << (int)cellptr->neighbors.size() << " " << (int)cellptr->gapneighbors.size() << " "
               << "\n";

    cellptr->gapneighbors.insert(gapneighborptr->globalid);
    gapneighborptr->gapneighbors.insert(cellptr->globalid);
    iadded++;
  }
}

//This function iterates over the edges of the cells and returns the minimum distance
double GeoExtractor::cellsDelta(Cell *cp1, Cell *cp2)
{
  DetId id1 = cp1->globalid;
  DetId id2 = cp2->globalid;
  LOG(DEBUG) << "Start GeoExtractor::cellsDelta for " << id1.rawId() << "\n";
  LOG(DEBUG) << *getCellPtr(id1) << "\n";
  LOG(DEBUG) << getCellHash(id1) << "\n";
  LOG(DEBUG) << "and " << id2.rawId() << "\n";
  LOG(DEBUG) << *getCellPtr(id2) << "\n";
  LOG(DEBUG) << getCellHash(id2) << "\n";

  edm::ESHandle<HGCalGeometry> &geo1 = m_geom[id1.det()];
  edm::ESHandle<HGCalGeometry> &geo2 = m_geom[id2.det()];

  LOG(DEBUG) << "id1 is present in geo1 " << geo1->present(id1) << "\n";
  LOG(DEBUG) << "id2 is present in geo2 " << geo2->present(id2) << "\n";

  LOG(DEBUG) << "area cell 1" << geo1->getArea(id1) << "\n";
  LOG(DEBUG) << "area cell 2" << geo2->getArea(id2) << "\n";

  std::vector<GlobalPoint> v_corners1 = geo1->getCorners(id1);
  std::vector<GlobalPoint> v_corners2 = geo2->getCorners(id2);
  GlobalPoint p1 = geo1->getPosition(id1);
  GlobalPoint p2 = geo2->getPosition(id2);

  LOG(DEBUG) << "id1 pos" << p1 << "\n";
  LOG(DEBUG) << "id2 pos" << p2 << "\n";

  LOG(DEBUG) << "centers diff: " << dist(p1, p2) << "\n";

  double foodist= 999999;
  double curdist;
  std::vector<GlobalPoint>::iterator bestc1;
  std::vector<GlobalPoint>::iterator bestc2;
  for (std::vector<GlobalPoint>::iterator c1it = v_corners1.begin(); c1it != v_corners1.end(); ++c1it) {
    for (std::vector<GlobalPoint>::iterator c2it = v_corners2.begin(); c2it != v_corners2.end(); ++c2it){
      curdist =dist(*c1it,*c2it);
      if (curdist<foodist){
        foodist = curdist;
        bestc1=c1it;
        bestc2=c2it;
      }
    }
  }
  LOG(DEBUG) << "MinimumDelta "<< foodist << " c1 " << (bestc1 - v_corners1.begin()) << *bestc1 << " c2 " << (bestc2 - v_corners2.begin()) << *bestc2 << "\n";

  return foodist;
  //The following algorithm solves the problem much faster then the previous double loop, but doesnt yield the minumum distance for irregular geometry, eg. around 2752652289.


  std::vector<GlobalPoint>::iterator cell1corner = v_corners1.begin();
  std::vector<GlobalPoint>::iterator cell2corner = v_corners2.begin();

  //Start witht the cells closest to  the other center:
  double olddist= dist(*cell1corner, p2);
  for (std::vector<GlobalPoint>::iterator i = v_corners1.begin()+1; i != v_corners1.end(); ++i)
    {
      if (dist(*i, p2) < olddist)
      cell1corner=i;
    }
  LOG(DEBUG)<< "cell 1 Center: "<< p1 << "Start edge for cell 1 " << *cell1corner << " at distance to the other center "<< olddist <<"\n" ;

  olddist= dist(*cell2corner, p1);
  for (std::vector<GlobalPoint>::iterator i = v_corners2.begin()+1; i != v_corners2.end(); ++i)
  {
    if (dist(*i, p1) < olddist)
    cell2corner=i;
  }
  LOG(DEBUG)<< "cell 2 Center: "<< p2 << "Start edge for cell 2 " << *cell2corner << " at distance to the other center "<< olddist <<"\n" ;

  //Find the minimum distance between the corners
  LOG(DEBUG) << "pre dist\n";
  double delta = dist(*cell1corner, *cell2corner);
  double tmpdelta;
  bool improvement = true;
  LOG(DEBUG) << "pre loop\n";
  while (improvement)
  {
    LOG(DEBUG) << "Corner1:" << *cell1corner << " Corner2: " << *cell2corner << " delta "<< delta <<"\n";
    improvement = false;
    //Corner1 forward
    if (cell1corner + 1 != v_corners1.end())
    {
      tmpdelta = dist(*(cell1corner + 1), *cell2corner);
      if (tmpdelta < delta)
      {
        cell1corner++;
        delta = tmpdelta;
        improvement = true;
        LOG(DEBUG) << "Corner1 ++ noloop c1 " << (cell1corner - v_corners1.begin()) << " c2 " << (cell2corner - v_corners2.begin()) << "\n";
        continue;
      }
    }
    //Corner1 vector loop
    else
    {
      tmpdelta = dist(*v_corners1.begin(), *cell2corner);
      if (tmpdelta < delta)
      {
        cell1corner = v_corners1.begin();
        delta = tmpdelta;
        improvement = true;
        LOG(DEBUG) << "Corner1 ++ loop c1 " << (cell1corner - v_corners1.begin()) << " c2 " << (cell2corner - v_corners2.begin()) << "\n";
        continue;
      }
    }
    //Corner1 backwards
    if (cell1corner != v_corners1.begin())
    {
      tmpdelta = dist(*(cell1corner - 1), *cell2corner);
      if (tmpdelta < delta)
      {
        cell1corner--;
        delta = tmpdelta;
        improvement = true;
        LOG(DEBUG) << "Corner1 -- noloop c1 " << (cell1corner - v_corners1.begin()) << " c2 " << (cell2corner - v_corners2.begin()) << "\n";
        continue;
      }
    }
    //Corner1 vector loop
    else
    {
      tmpdelta = dist(*(v_corners1.end() - 1), *cell2corner);
      if (tmpdelta < delta)
      {
        cell1corner = (v_corners1.end() - 1);
        delta = tmpdelta;
        improvement = true;
        LOG(DEBUG) << "Corner1 -- loop c1 " << (cell1corner - v_corners1.begin()) << " c2 " << (cell2corner - v_corners2.begin()) << "\n";
        continue;
      }
    }
    //
    //Conver vector 2
    //
    if (cell2corner + 1 != v_corners2.end())
    {
      tmpdelta = dist(*(cell2corner + 1), *cell1corner);
      if (tmpdelta < delta)
      {
        cell2corner++;
        delta = tmpdelta;
        improvement = true;
        LOG(DEBUG) << "Corner2 ++ noloop c1 " << (cell1corner - v_corners1.begin()) << " c2 " << (cell2corner - v_corners2.begin()) << "\n";
        continue;
      }
    }
    //Corner1 vector loop
    else
    {
      tmpdelta = dist(*v_corners2.begin(), *cell1corner);
      if (tmpdelta < delta)
      {
        cell2corner = v_corners2.begin();
        delta = tmpdelta;
        improvement = true;
        LOG(DEBUG) << "Corner2 ++ loop c1 " << (cell1corner - v_corners1.begin()) << " c2 " << (cell2corner - v_corners2.begin()) << "\n";
        continue;
      }
    }
    //Corner1 backwards
    if (cell2corner != v_corners2.begin())
    {
      tmpdelta = dist(*(cell2corner - 1), *cell1corner);
      if (tmpdelta < delta)
      {
        cell2corner--;
        delta = tmpdelta;
        improvement = true;
        LOG(DEBUG) << "Corner2 -- noloop c1 " << (cell1corner - v_corners1.begin()) << " c2 " << (cell2corner - v_corners2.begin()) << "\n";
        continue;
      }
    }
    //Corner1 vector loop
    else
    {
      tmpdelta = dist(*(v_corners2.end() - 1), *cell1corner);
      if (tmpdelta < delta)
      {
        cell2corner = (v_corners2.end() - 1);
        delta = tmpdelta;
        improvement = true;
        LOG(DEBUG) << "Corner2 -- loop c1 " << (cell1corner - v_corners1.begin()) << " c2 " << (cell2corner - v_corners2.begin()) << "\n";
        continue;
      }
    }
  }
  LOG(DEBUG) << "final delta " << delta << "\n\n";
  return delta;
}

bool GeoExtractor::rangecond(std::vector<PosListTup>::iterator iter, Cell *cellptr)
{
  return (std::abs(std::get<0>(*iter) - cellptr->x) < maxSearchDelta);
}

// double GeoExtractor::xposdiff(std::vector<PosListTup>::iterator iter, Cell *cellptr)
// {
//   return (std::abs(std::get<0>(*iter) - cellptr->x));
// }
double GeoExtractor::xposdiff(std::vector<PosListTup>::iterator iter, Cell *cellptr)
{
  return (std::abs(std::get<1>(*iter)->x - cellptr->x));
}

double GeoExtractor::dist(GlobalPoint &p1, GlobalPoint &p2)
{
  double dx = (p1.x() - p2.x());
  double dy = (p1.y() - p2.y());
  double dz = (p1.z() - p2.z());
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}
