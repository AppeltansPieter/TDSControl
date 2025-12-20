classdef TDS < handle

    properties
        m_TDS
    end

    methods
        function obj = TDS(A, hA)
           obj.m_TDS = wrap_tds_create(A, hA);
        end

        function l = tds_roots(obj, N)
            l = wrap_tds_roots(obj.m_TDS, N);
        end
    end
end
