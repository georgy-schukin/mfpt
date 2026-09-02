#pragma once

#include "didal/data/sync_distributed_storage.h"
#include "didal/base/index.h"

#include "shadowed_array2d.h"

template <typename T>
class DistributedMesh2D : public ddl::SyncDistributedStorage<ddl::Index<2>, ShadowedArray2D<T>> {
public:
private:
};
