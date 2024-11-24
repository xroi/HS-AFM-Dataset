import copy


class FGParams:
    def __init__(self,
                 res_from=None,
                 res_to=None,
                 self_k=None,
                 self_range=None,
                 kap_k=None,
                 kap_range=None,
                 nonspec_k=None,
                 nonspec_range=None,
                 backbone_k=None,
                 backbone_tau=None):
        #    assert(res_to>=res_from)
        self.res_from = res_from
        self.res_to = res_to
        self.self_k = self_k
        self.self_range = self_range
        self.kap_k = kap_k
        self.kap_range = kap_range
        self.nonspec_k = nonspec_k
        self.nonspec_range = nonspec_range
        self.backbone_k = backbone_k
        self.backbone_tau = backbone_tau

    def deepcopy(self):
        return copy.deepcopy(self)

    def get_copy(self, res_from=None, res_to=None, self_k=None, nonspec_k=None):
        ret = self.deepcopy()
        ret.res_from = res_from
        ret.res_to = res_to
        if self_k is not None:
            ret.self_k = self_k
        if nonspec_k is not None:
            ret.nonspec_k = nonspec_k
        return ret

    def update(self, dictionary):
        for key, value in dictionary.items():
            setattr(self, key, value)

    def __str__(self):
        ''' Returns a string representation of an FGParams object '''
        return f"FGParams(res_from={self.res_from}, res_to={self.res_to}," + \
            f"\n\tself_k={self.self_k}, self_range={self.self_range},\n\tkap_k={self.kap_k}," + \
            f" kap_range={self.kap_range},\n\tnonspec_k={self.nonspec_k}," + \
            f" nonspec_range={self.nonspec_range},\n\tbackbone_k={self.backbone_k}," + \
            f" backbone_tau={self.backbone_tau})"