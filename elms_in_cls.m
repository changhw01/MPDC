function I_res = elms_in_cls(cls,I_version)

n_versions = length(I_version);

I_res = [];
for i = 1:n_versions
    I_new = find(cls.rep_clst{I_version(i)});
    I_res = cat(2,I_res,setdiff(I_new,I_res));
end
