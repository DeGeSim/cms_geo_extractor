#pragma once

std::vector<DetId> GeoExtractor::filterCellIds(const std::vector<DetId> v_allCellIds)
{
  std::vector<DetId> v_validHGCalIds;
  LOG(DEBUG) << "#All cell ids " << v_allCellIds.size() << "\n";
  for (auto detectorid : v_HGCalDets)
  {
    m_rej[detectorid];
    m_rej[detectorid]["det"] = 0;
    m_rej[detectorid]["x"] = 0;
    m_rej[detectorid]["y"] = 0;
    m_rej[detectorid]["z"] = 0;
  }
  for (int i = 0; i < (int)v_allCellIds.size(); i++)
  {
    auto detectorid = v_allCellIds[i].det();
    // Skip IDs from other detector parts
    if (detectorid != DetId::HGCalEE && detectorid != DetId::HGCalHSi &&
        detectorid != DetId::HGCalHSc)
    {
      continue;
    }
    m_rej[detectorid]["det"]++;
    // Radius = 25

    auto x = recHitTools.getPosition(v_allCellIds[i]).x();
    if (x < -50 || x > 50)
    {
      continue;
    }
    m_rej[detectorid]["x"]++;

    auto y = recHitTools.getPosition(v_allCellIds[i]).y();
    if (y < -150 || y > -50)
    {
      continue;
    }
    m_rej[detectorid]["y"]++;

    auto z = recHitTools.getPosition(v_allCellIds[i]).z();
    if (z < 0)
    {
      continue;
    }
    m_rej[detectorid]["z"]++;
    // If all checks are pass, add the detId to the list of valid Ids.
    v_validHGCalIds.push_back(v_allCellIds[i]);
  }
  LOG(INFO) << "\tdet    ";
  LOG(INFO) << "x    ";
  LOG(INFO) << "y    ";
  LOG(INFO) << "z \n";
  for (auto detectorid : v_HGCalDets)
  {
    LOG(INFO) << detectorid << "\t";
    LOG(INFO) << m_rej[detectorid]["det"] << "\t";
    LOG(INFO) << m_rej[detectorid]["x"] << "\t";
    LOG(INFO) << m_rej[detectorid]["y"] << "\t";
    LOG(INFO) << m_rej[detectorid]["z"] << "\n";
  }
  LOG(INFO) << "Cells left: " << v_validHGCalIds.size() << "\n";
  return v_validHGCalIds;
}



bool GeoExtractor::validId(DetId id)
{
  auto [detectorid, subdetid, layerid, waferortileid, cellid] = getCellHash(id);
  // unsigned int detectorid = id.det();
  // unsigned int subdetid = id.subdetId();
  // unsigned int layerid = recHitTools.getLayer(id);
  if (detectorid != DetId::HGCalEE && detectorid != DetId::HGCalHSi && detectorid != DetId::HGCalHSc)
  {
    LOG(DEBUG) << "wrong detector " << detectorid << "\n";
    return false;
  }
  if (subdetid != 0)
  {
    LOG(DEBUG) << "wrong subdetector " << subdetid << "\n";
    return false;
  }
  // if (detectorid != DetId::HGCalEE)
  // {
  //   if (layerid < 1 || layerid > 28)
  //   {
  //     LOG(DEBUG) << "layer ee "<< layerid << "\n";
  //     return false;
  //   }
  // }
  // if (detectorid != DetId::HGCalHSi)
  // {
  //   if (layerid < 1 || layerid > 22)
  //   {
  //     LOG(DEBUG) << "layer si "<< layerid << "\n";
  //     return false;
  //   }
  // }
  // if (detectorid != DetId::HGCalHSc)
  // {
  //   if (layerid < 9 || layerid > 22)
  //   {
  //     LOG(DEBUG) << "layer hsc "<< layerid << "\n";
  //     return false;
  //   }
  // }
  return true;
}
void GeoExtractor::validateId(DetId id)
{
  if (!validId(id))
  {
    LOG(DEBUG) << "Invalid cell Id:" << id.rawId()<<"\n";
    LOG(DEBUG) << *getCellPtr(id)<<"\n";
    throw std::invalid_argument("Invalid Cell");
  }
}
bool GeoExtractor::isSiliconDet(int ndet)
{
  if (ndet == DetId::HGCalEE || ndet == DetId::HGCalHSi)
  {
    return true;
  }
  else if (ndet == DetId::HGCalHSc)
  {
    return false;
  }
  else
  {
    throw std::invalid_argument("Invalid detector");
    return false;
  }
}