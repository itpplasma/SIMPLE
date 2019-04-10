function neo_orb_init(neo_orb_root)
    addpath([neo_orb_root, '/lib'])

    if not(libisloaded('libneo_orb'))
        loadlibrary('libneo_orb', [neo_orb_root, '/matlab/neo_orb.h'])
    end

    libfunctions('libneo_orb', '-full')
end
