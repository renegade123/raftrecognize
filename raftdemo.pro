pro raftdemo
obj = obj_new("raftrecognize",'C:\Users\name\IDLWorkspace83\raftrecognize\data\testme.bmp',1)
obj->superclassfy
obj->filterPlaque
obj->dilateerode
obj->raster_to_vector
end