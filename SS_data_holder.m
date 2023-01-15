classdef SS_data_holder < handle & matlab.mixin.Heterogeneous
    
    % ============ Fields =================================
    
    properties
        
        mouse_name
                
        freqs % list of all presented sound frequencies (kHz)
        bandwidths % list of all presented bandwidths (% oct.)
        dur % sound duration (ms)
        n_trials % number of presentations of wach stiimulus
        ISI % mean inter stimulus intervals (seconds)
        ISI_range % range of inter stimulus interval around the mean
        seq % sequence of presentation of stimuli
                
        area_names % list of brain areas of each neurons
        area_acronyms % acronyms of brain areas
        sorted_rasters % raster plots for BBS
        sorted_rasters2 % raster plots for notches
        PSTH_bands % PSTHs for all BBS
        PSTH_notches % PSTH for all notches
        time_vec % time vector for plotting
        sig_cells_bands % indicator of significance responses to various BBS response phases
        sig_cells_notches % indicator of significance responses to various notch response phases
        amps_bands % response amplitudes for BBS
        amps_notches % response amplitudes for notches
        shuff_amps_bands % shuffled response amplitudes for BBS
        
        bands_mat % matrix of BBS
        notches_mat % matris of notches
        n_exp_freqs % number of frequencies oin experiment
        n_bands % number of bandwidths in experiment
        n_sim_freqs % number of frequencies in model
        sim_freqs % frequencies in model
        total_oct % number of octaves in model
                        
        fit_parameters % fit parameters
        fit_parameters_for_disp % fit parametersfor display
        fit_corrs % BBS fit correlations
        fit_corrs_shuff % BBS fit correlations for shuffled responses
        notch_corrs % notch fit correlations 
        notch_corrs_shuff % notch fit correlations for shuffled responses
        
        top_freq_fit % top frequency gaussian fit parameters
        top_freq_R % top frequency gaussian fit R^2
        
        bottom_freq_fit % bottom frequency gaussian fit parameters
        bottom_freq_R % bottom frequency gaussian fit R^2
        
    end
    
    properties (Constant = true)
        
        BIN_SIZE = 10 ; % bin size for PSTH
        PRE_STIM_TIME = 100 ; % ms before stim onset
        POST_STIM_TIME = 100 ; % ms after stim onset
                
        oct_padding = 3 ; % 
        
        filter_power = 1 ; % shape of cochlear filter
        filter_width = 50 ; % width of cochlear filter
        
        LM = [0 0.05 0 0] ; % lagrange multipliers for fit (only used to limit width)
        
        n_shuffle_itter = 100 ; % number of otteration for shuffle
                
    end
    
    % ============ Methods ================================
    
    methods
        
        function obj = SS_data_holder(mouse_name)
            
            obj.mouse_name = mouse_name ;
            
        end 
              
        function organize_data(obj)
            
            obj.n_exp_freqs = length(obj.freqs) ;
            obj.n_bands = length(obj.bandwidths) ;
            obj.n_sim_freqs = 100*(obj.n_exp_freqs+8*obj.oct_padding-1)+1 ;
            obj.sim_freqs = logspace(log10(obj.freqs(1)/(2^obj.oct_padding)),log10(obj.freqs(end)*2^obj.oct_padding),obj.n_sim_freqs) ;
            obj.total_oct = log2(obj.freqs(end)/obj.freqs(1)) + 2*obj.oct_padding ;
            
%             obj.make_PSTH('bands') 
%             obj.check_sig_cells('bands')
%             obj.make_bands_matrix()
%             obj.fit_neurons()
%             obj.fit_neurons('shuffle')
%             obj.top_freq_tuning()
%             obj.bottom_freq_tuning()
%             
%             obj.make_notches_matrix()
%             obj.make_PSTH('notches')
%             obj.check_sig_cells('notches')
%             obj.predict_notches()

                           
        end
                
        function make_PSTH(obj,type)
                      
            PSTH = {} ;
            switch type
                case 'bands'
                    sorted_rasters = obj.sorted_rasters ;
                case 'notches'
                    sorted_rasters = obj.sorted_rasters2 ;
            end
                        
            for neuron = 1 : size(sorted_rasters,1) 
                for freq = 1:size(sorted_rasters,2)
                    for bandwidth = 1:size(sorted_rasters,3)
                        rasters = [] ;
                        for trial = 1 : size(sorted_rasters{neuron,freq,bandwidth},1)
                            temp1 = sorted_rasters{neuron,freq,bandwidth}(trial,:) ;
                            temp2 = sum(reshape(temp1,obj.BIN_SIZE,length(temp1)/obj.BIN_SIZE),1) ;
                            rasters = [rasters ; temp2] ;
                        end
                        PSTH{neuron,freq,bandwidth} = mean(rasters,1)*1000/obj.BIN_SIZE ;
                    end
                end
            end
            
            switch type
                case 'bands'
                    obj.PSTH_bands = PSTH ;
                    temp = length(obj.PSTH_bands{1,1,1}) ;
                    obj.time_vec = obj.BIN_SIZE*(0:(temp-1)) - obj.PRE_STIM_TIME ;
                case 'notches'
                    obj.PSTH_notches = PSTH ;
            end
                       
        end
        
        function check_sig_cells(obj,type)
            
            switch type
                case 'bands' 
                    sorted_rasters = obj.sorted_rasters ;
                case 'notches'
                    sorted_rasters = obj.sorted_rasters2 ;
            end
            
            sig_stims = zeros(size(sorted_rasters)) ;
            sig_cells = zeros(4, size(sorted_rasters,1)) ;
            amps = zeros(size(sorted_rasters)) ;
            
            pre_window = obj.PRE_STIM_TIME/obj.BIN_SIZE ;
            post_window = obj.POST_STIM_TIME/obj.BIN_SIZE ;
            
            a = 0.05 / (obj.n_exp_freqs * (obj.n_bands+1)) ;
            
            for neuron = 1 : size(sorted_rasters,1)
                amps_all = [] ;
                off_amps_all = [] ;
                for freq = 1:size(sorted_rasters,2)
                    for bandwidth = 1:size(sorted_rasters,3)

                        rasters = [] ;
                        for trial = 1 : size(sorted_rasters{neuron,freq,bandwidth},1)
                            temp1 = sorted_rasters{neuron,freq,bandwidth}(trial,:) ;
                            temp2 = sum(reshape(temp1,obj.BIN_SIZE,length(temp1)/obj.BIN_SIZE),1) ;
                            rasters = [rasters ; temp2] ;
                        end

                        pre_stim = sum(rasters(:,(11-pre_window):10),2)/(obj.PRE_STIM_TIME/1000) ;
                        post_stim = sum(rasters(:,10 + (1:post_window)),2)/(obj.POST_STIM_TIME/1000) ;
                        off_stim = sum(rasters(:,10 + post_window +  (1:post_window)),2)/(obj.POST_STIM_TIME/1000) ;


                        amp = mean(post_stim - pre_stim) ;
                        amps(neuron,freq,bandwidth) = amp ;
                        amps_all = [amps_all ; post_stim - pre_stim] ;
                        off_amps_all = [off_amps_all ; off_stim - pre_stim] ;

                    end
                end
                sig = [ttest(amps_all,0,'Tail','Right','Alpha',a), ...
                       ttest(amps_all,0,'Tail','Left','Alpha',a), ...
                       ttest(off_amps_all,0,'Tail','Right','Alpha',a), ...
                       ttest(off_amps_all,0,'Tail','Left','Alpha',a)];
                sig(isnan(sig)) = 0 ;
                sig_cells(:,neuron) = sig ;
            end
            
            switch type
                case 'bands' 
                    obj.sig_cells_bands = sig_cells ;
                    obj.amps_bands = amps ;
                    obj.shuff_amps_bands = amps ;
                case 'notches'
                    obj.sig_cells_notches = sig_cells ;
                    obj.amps_notches = amps ;
            end
            
        end
                
        function make_bands_matrix(obj)
               
            bands = [0 obj.bandwidths] ;
            obj.bands_mat = zeros(obj.n_exp_freqs * (obj.n_bands + 1) , obj.n_sim_freqs) ;
            
            for i = 1 : obj.n_exp_freqs
                for j = 1 : length(bands)
                    freq = obj.freqs(i) ;
                    band = bands(j) ;
                    temp = zeros(1,obj.n_sim_freqs) ;
                    temp(obj.sim_freqs > (freq*2^(-band/2))*0.9999 & obj.sim_freqs < (freq*2^(band/2))*1.0001) = 1 ;
                    temp = obj.cochlear_filter(temp) ;
                    obj.bands_mat((i-1)*(obj.n_bands + 1) + j , :) = temp ;
                end
            end
            
        end
        
        function make_notches_matrix(obj)
               
            bands = obj.bandwidths ;
            notches_mat1 = zeros(obj.n_exp_freqs * obj.n_bands , obj.n_sim_freqs) ;
            notches_mat2 = zeros(obj.n_exp_freqs * obj.n_bands , obj.n_sim_freqs) ;
            
            for i = 1 : obj.n_exp_freqs
                for j = 1 : length(bands)
                    freq = obj.freqs(i) ;
                    band = bands(j) ;
                    temp1 = zeros(1,obj.n_sim_freqs) ;
                    temp1(obj.sim_freqs < (freq*2^(-band/2))*0.9999) = 1 ;
                    temp2 = zeros(1,obj.n_sim_freqs) ;
                    temp2(obj.sim_freqs > (freq*2^(band/2))*1.0001) = 1 ;
                    temp1 = obj.cochlear_filter(temp1 , 'left to notch') ;
                    temp2 = obj.cochlear_filter(temp2 , 'right to notch') ;
                   notches_mat1((i-1)*obj.n_bands + j , :) = temp1 ;
                    notches_mat2((i-1)*obj.n_bands + j , :) = temp2 ;
                end
            end
            
            obj.notches_mat = max(cat(3,notches_mat1,notches_mat2),[],3) ;
            
        end
        
        function tuning = make_tuning_curve(obj,Mean,Width,Skew,IE)

            x = 1 : obj.n_sim_freqs ;
            
            y1 = (Width^2 -IE*(x - Mean).^2) .* exp( ( -(x - Mean).^2 ) / (2*Width^2) ) ;
            y2 = Skew*normcdf(x,Mean,Width/2) + (1-Skew)*(1-normcdf(x,Mean,Width/2)) ;

            tuning = y1.*y2 ;
            tuning = tuning/max(tuning) ;

        end
        
        function responses = calc_sim_band_responses(obj,Mean,Width,Skew,IE)
                        
            tuning = obj.make_tuning_curve(Mean,Width,Skew,IE) ;
            
            temp = rectify(obj.bands_mat*tuning') ;

            responses = reshape(temp,obj.n_bands + 1,obj.n_exp_freqs)' ;
            
        end
        
        function responses = calc_sim_notch_responses(obj,Mean,Width,Skew,IE)
            
            
            tuning = obj.make_tuning_curve(Mean,Width,Skew,IE) ;
            
            temp = rectify(obj.notches_mat*tuning') ;

            responses = reshape(temp,obj.n_bands,obj.n_exp_freqs)' ;
            
        end
        
        function f_sig = cochlear_filter(obj,sig,side)
            
            if nargin < 3
                side = 'band' ;
            end
            
            n_points = obj.n_sim_freqs * (obj.filter_width/100) / obj.total_oct ;
            
            temp1 = linspace(0,1,floor(1*n_points)).^obj.filter_power ;
            temp2 = fliplr(linspace(0,1,ceil(1*n_points))).^obj.filter_power ;
            a = find(sig,1,'first') ;
            b = find(sig,1,'last') ;
            
            f_sig = sig ;
            
            if ~strcmp(side,'left to notch')
                f_sig((a-length(temp1)) : (a-1)) = temp1 ; 
            end
            
            if ~strcmp(side,'right to notch')
                f_sig((b+1) : (b+length(temp2))) = temp2 ;
            end
                 
        end
        
        function corr_val = calc_mexican(obj,neuron,x,fit_mode)
            
            if nargin == 4 && strcmp(fit_mode,'shuffle')
                r = squeeze(obj.shuff_amps_bands(neuron,:,2:end)) ;
            else
                r = squeeze(obj.amps_bands(neuron,:,2:end)) ;
            end
            Mean = x(1) ; Width = x(2) ; Skew = x(3) ; IE = x(4) ;
            r_hat = obj.calc_sim_band_responses(Mean,Width,Skew,IE) ;
            temp = corrcoef(r,r_hat(:,2:end)) ;
            x_hat(1) = obj.sim_freqs(floor(x(1))) ;
            x_hat(2) = obj.total_oct*x(2)/obj.n_sim_freqs ;
            x_hat(3) = -1+2 * x(3) ;
            x_hat(4) = x(4) ;
            corr_val = obj.LM*abs(x_hat)'-temp(1,2) ;
    
        end
        
        function x_hat = fit_mexican(obj,neuron,fit_mode)
            
            if nargin == 3 && strcmp(fit_mode,'shuffle')
                r = squeeze(obj.shuff_amps_bands(neuron,:,2:end)) ;
                fun = @(x)obj.calc_mexican(neuron,x,fit_mode) ;
            else
                r = squeeze(obj.amps_bands(neuron,:,2:end)) ;
                fun = @(x)obj.calc_mexican(neuron,x) ;
            end
            
            Mean_min = find(obj.sim_freqs > obj.freqs(1) , 1 , 'first') ;
            Width_min = 0.02*obj.n_sim_freqs/obj.total_oct ;
            Skew_min = 0 ;
            IE_min = 0 ;

            Mean_max = find(obj.sim_freqs < obj.freqs(end) , 1 , 'last') ;
            Width_max = obj.total_oct*obj.n_sim_freqs/obj.total_oct ;
            Skew_max = 1 ;
            IE_max = 1 ;

            [~,temp] = max(mean(r,2)) ;
            Mean_init = find(obj.sim_freqs > obj.freqs(temp), 1 , 'first') ;
            Width_init = (Width_min + Width_max)/2 ;
            Skew_init = (Skew_min + Skew_max)/2 ;
            IE_init = (IE_min + IE_max)/2 ;

            options = optimset('MaxFunEvals',3000,'MaxIter',3000) ;
            x_hat = fminsearchbnd(fun,[Mean_init Width_init Skew_init IE_init],...
                [Mean_min Width_min Skew_min IE_min],[Mean_max Width_max Skew_max IE_max],options) ;


        end

        function fit_neurons(obj,fit_mode)
            
            for neuron = 1 : size(obj.sorted_rasters,1)
                if obj.sig_cells_bands(1,neuron)
                    if nargin == 2 && strcmp(fit_mode,'shuffle')
                        for i = 1 : obj.n_shuffle_itter
                            temp = obj.shuff_amps_bands(neuron,:,2:end) ;
                            temp = temp(:) ;
                            temp = temp(randperm(numel(temp))) ;
                            obj.shuff_amps_bands(neuron,:,2:end) = reshape(temp,size(obj.shuff_amps_bands,2),size(obj.shuff_amps_bands,3)-1) ;
                            temp = fit_mexican(obj,neuron,fit_mode) ;
                            x_hat(1) = obj.sim_freqs(floor(temp(1))) ;
                            x_hat(2) = obj.total_oct*temp(2)/obj.n_sim_freqs ;
                            x_hat(3) = -1+2 * temp(3) ;
                            x_hat(4) = temp(4) ;
                            Corrs(i) = obj.LM*abs(x_hat)' - calc_mexican(obj,neuron,temp,fit_mode) ;
                        end
                        obj.fit_corrs_shuff(neuron) = mean(Corrs) ;
                    else
                        temp = fit_mexican(obj,neuron) ;
                        x_hat(1) = obj.sim_freqs(floor(temp(1))) ;
                        x_hat(2) = obj.total_oct*temp(2)/obj.n_sim_freqs ;
                        x_hat(3) = -1+2 * temp(3) ;
                        x_hat(4) = temp(4) ;
                        obj.fit_parameters(neuron,:) = x_hat ;
                        obj.fit_parameters_for_disp(neuron,:) = temp ;
                        obj.fit_corrs(neuron) = obj.LM*abs(x_hat)' - calc_mexican(obj,neuron,temp) ;
                    end
                else
                    obj.fit_parameters(neuron,:) = zeros(1,4) ;
                    obj.fit_parameters_for_disp(neuron,:) = zeros(1,4) ;
                    obj.fit_corrs(neuron) = 0 ;
                    obj.fit_corrs_shuff(neuron) = 0 ;
                end
            end
            
        end
        
        function predict_notches(obj)
            
            obj.notch_corrs = zeros(1,size(obj.sorted_rasters,1)) ;
            
            for neuron = 1 : size(obj.sorted_rasters,1)
                if obj.sig_cells_bands(1,neuron)
                    r = squeeze(obj.amps_notches(neuron,:,2:end))' ;
                    
                    x_hat = obj.fit_parameters_for_disp(neuron,:) ;    
                    r_hat = obj.calc_sim_notch_responses(x_hat(1),x_hat(2),x_hat(3),x_hat(4))' ;
                    temp = corrcoef(r,r_hat) ;
                    obj.notch_corrs(neuron) = temp(1,2) ;
                    shuff_corrs = [] ;
                    for i = 1 : obj.n_shuffle_itter
                        temp = r(:) ;
                        temp = temp(randperm(numel(temp))) ;
                        r_shuff = reshape(temp,size(r,1),size(r,2)) ;
                        temp = corrcoef(r_shuff,r_hat) ;
                        shuff_corrs(i) = temp(1,2) ;          
                    end
                    obj.notch_corrs_shuff(neuron) = mean(shuff_corrs) ;
                end
            end
            
            
        end
        
        function top_freq_tuning(obj)
            
            obj.top_freq_R = [] ;
            
            for neuron = 1 : size(obj.sorted_rasters,1)
                if obj.sig_cells_bands(1,neuron)
                    resps = squeeze(obj.amps_bands(neuron,:,2:end)) ;
                    top_freqs = zeros(length(obj.freqs),length(obj.bandwidths)) ;
                    for f = 1 : length(obj.freqs)
                        for b = 1 : length(obj.bandwidths)
                            top_freqs(f,b) = obj.freqs(f)*2^(obj.bandwidths(b)/2) ;
                        end
                    end
                    x = top_freqs(:) ;
                    y = resps(:) ;
                    try
                        [f,gof] = fit(log(x),y,'gauss1') ;
                        obj.top_freq_fit(neuron,:) = [f.a1 f.b1 f.c1] ;
                        obj.top_freq_R(neuron) = gof.rsquare ;
                    catch
                        obj.top_freq_fit(neuron,:) = zeros(1,3) ;
                        obj.top_freq_R(neuron) = 0 ;
                    end
                else
                    obj.top_freq_fit(neuron,:) = zeros(1,3) ;
                    obj.top_freq_R(neuron) = 0 ;
                end
            end
            
        end
        
        function bottom_freq_tuning(obj)
            
            obj.bottom_freq_R = [] ;
            
            for neuron = 1 : size(obj.sorted_rasters,1)
                if obj.sig_cells_bands(1,neuron)
                    resps = squeeze(obj.amps_bands(neuron,:,2:end)) ;
                    bottom_freqs = zeros(length(obj.freqs),length(obj.bandwidths)) ;
                    for f = 1 : length(obj.freqs)
                        for b = 1 : length(obj.bandwidths)
                            bottom_freqs(f,b) = obj.freqs(f)*2^(-obj.bandwidths(b)/2) ;
                        end
                    end
                    x = bottom_freqs(:) ;
                    y = resps(:) ;
                    try
                        [f,gof] = fit(log(x),y,'gauss1') ;
                        obj.bottom_freq_fit(neuron,:) = [f.a1 f.b1 f.c1] ;
                        obj.bottom_freq_R(neuron) = gof.rsquare ;
                    catch
                        obj.bottom_freq_fit(neuron,:) = zeros(1,3) ;
                        obj.bottom_freq_R(neuron) = 0 ;
                    end
                        
                else
                    obj.bottom_freq_fit(neuron,:) = zeros(1,3) ;
                    obj.bottom_freq_R(neuron) = 0 ;
                end
            end
            
        end
                       
    end
    
end