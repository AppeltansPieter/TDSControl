module TDSControl
    using CxxWrap
    @wrapmodule(() -> joinpath("/home/pietera/workspace/tds-control-v2/build", "libtdscontroljl.so"))


    function __init__()
        @initcxx
    end
end