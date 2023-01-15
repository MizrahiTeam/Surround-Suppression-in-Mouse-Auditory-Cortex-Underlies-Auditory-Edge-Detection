% ========================== Opening / closing functions =================

function varargout = SS_anlss(varargin)

    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @SS_anlss_OpeningFcn, ...
                       'gui_OutputFcn',  @SS_anlss_OutputFcn, ...
                       'gui_LayoutFcn',  [], ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
       gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end

end

function SS_anlss_OpeningFcn(hObject, eventdata, handles, varargin)

    handles.output = hObject;
    handles.min_corr = 0.644 ;
    set(handles.figure1,'windowscrollWheelFcn',@scroll_func)

    guidata(hObject, handles);
    
end

function varargout = SS_anlss_OutputFcn(hObject, eventdata, handles)

    varargout{1} = handles.output;
    
end

% ========================== Main functions ==============================

function anlss_type_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    Ax = handles.anlss_ax ;
    cla(handles.anlss_ax , 'reset')
    axis(handles.anlss_ax , 'auto')
    hold(handles.anlss_ax,'off')

    switch handles.anlss_type_menu.String{handles.anlss_type_menu.Value}
        case 'PSTH - bands'
            handles.anlss_type_menu.UserData = 1 ;
            plot_PSTH(Ax,hObject,'Bands')
        case 'PSTH - notches'
            handles.anlss_type_menu.UserData = 1 ;
            plot_PSTH(Ax,hObject,'Notches')
        case 'FBRA'
            handles.anlss_type_menu.UserData = 1 ;
            plot_FBRA(hObject,Ax)
        case 'nFBRA'
            handles.anlss_type_menu.UserData = 1 ;
            plot_nFBRA(hObject,Ax)
        case 'Tuning curve'
            handles.anlss_type_menu.UserData = 1 ;
            plot_tuning_curve(hObject,Ax)
        case 'Fitted FBRA'
            handles.anlss_type_menu.UserData = 1 ;
            plot_fitted_FBRA(hObject,Ax)
        case 'Predicted nFBRA'
            handles.anlss_type_menu.UserData = 1 ;
            plot_predicted_nFBRA(hObject,Ax)
        case 'Top frequency tuning'
            handles.anlss_type_menu.UserData = 1 ;
            plot_top_freq_tuning(hObject,Ax)
        case 'Bottom frequency tuning'
            handles.anlss_type_menu.UserData = 1 ;
            plot_bottom_freq_tuning(hObject,Ax)
        case 'Fit correlations dist.'
            handles.anlss_type_menu.UserData = 0 ;
            plot_fit_corrs_dist(hObject,Ax)
        case 'Notch correlations dist.'
            handles.anlss_type_menu.UserData = 0 ;
            plot_notch_corrs_dist(hObject,Ax)
        case 'Corrs scatter'
            handles.anlss_type_menu.UserData = 0 ;
            plot_corrs_scatter(hObject,Ax)
        case 'Mean values dist.'
            handles.anlss_type_menu.UserData = 0 ;
            plot_param_dist(hObject,Ax,1)
        case 'Width values dist.'
            handles.anlss_type_menu.UserData = 0 ;
            plot_param_dist(hObject,Ax,2)
        case 'Skew values dist.'
            handles.anlss_type_menu.UserData = 0 ;
            plot_param_dist(hObject,Ax,3)
        case 'I/E values dist.'
            handles.anlss_type_menu.UserData = 0 ;
            plot_param_dist(hObject,Ax,4)
        case 'Skew vs. I/E'
            handles.anlss_type_menu.UserData = 0 ;
            plot_skew_vs_IE(hObject,Ax)
        case 'Skew vs. I/E Q dist.'
            handles.anlss_type_menu.UserData = 0 ;
            plot_skew_vs_IE_quad_dist(hObject,Ax)
        case 'dR per quad'
            handles.anlss_type_menu.UserData = 0 ;
            plot_dR_fit_per_quad(hObject,Ax)
        case 'Top-freq vs. bottom-freq'
            handles.anlss_type_menu.UserData = 0 ;
            plot_top_freq_vs_bottom_freq(hObject,Ax)
            
    end
        
end

function update_anlss_type_menu(hObject)

    handles = guidata(hObject) ;
    handles.anlss_type_menu.String = {'PSTH - bands','PSTH - notches','FBRA','Fitted FBRA','nFBRA',...
        'Predicted nFBRA','Tuning curve','Top frequency tuning','Bottom frequency tuning',...
        'Fit correlations dist.','Notch correlations dist.',...
        'Corrs scatter','Mean values dist.','Width values dist.','Skew values dist.','I/E values dist.',...
        'Skew vs. I/E','Skew vs. I/E Q dist.','dR per quad','Top-freq vs. bottom-freq'} ;
    handles.anlss_type_menu.Value = 1 ;
    
end

function initialize_all(hObject)

    handles = guidata(hObject) ;
    
    handles.min_corr = 0.644 ;
    guidata(hObject,handles)
    for i = 1 : numel(handles.data_holder)
        handles.data_holder{i}.organize_data() ;
    end

    
    handles.anlss_type_menu.Enable = 'on' ;
    handles.mouse_menu.Enable = 'on' ;
    handles.neuron_menu.Enable = 'on' ;
    
    update_anlss_type_menu(hObject)
    mouse_menu_Callback(hObject,0,0)
    anlss_type_menu_Callback(hObject,0,0)
    
    guidata(hObject,handles)
        
end

% ========================== Plotting functions ==========================

function plot_PSTH(Ax,hObject,stim_type)

    handles = guidata(hObject) ;
    
    hold(Ax,'on')
        
    d = handles.data_holder{handles.mouse_menu.Value} ;
    switch stim_type
        case 'Bands'
            PSTH = d.PSTH_bands ;
        case 'Notches'
            PSTH = d.PSTH_notches ;
    end
    t = d.time_vec ;
    
    neuron = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;

    Ylims = [-0.1 0.1] ;
    
    for freq = 1:size(PSTH,2)
        for bandwidth = 1:size(PSTH,3)
            y{freq,bandwidth} = PSTH{neuron,freq,bandwidth} ;
            Ylims = [min([Ylims(1) min(y{freq,bandwidth})]) max([Ylims(2) max(y{freq,bandwidth})])] ;
        end
    end
    
    plot(Ax,[0 0],Ylims,'k')
    plot(Ax,[0 400],Ylims(1)*[1 1],'k')

    for freq = 1:size(PSTH,2)
        for bandwidth = 1:size(PSTH,3)
            plot(Ax,t + 1.5*freq*t(end),y{freq,bandwidth} - 1.1*(bandwidth-1)*diff(Ylims),'k')
        end
    end
    
    set(Ax,'XTick',[],'YTick',[])

end

function plot_top_freq_tuning(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder{handles.mouse_menu.Value} ;
    neuron = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;
    freqs = d.freqs ;
    bands = d.bandwidths ;
    resps = squeeze(d.amps_bands(neuron,:,2:end)) ;
    top_freqs = zeros(length(freqs),length(bands)) ;
    for f = 1 : length(freqs)
        for b = 1 : length(bands)
            top_freqs(f,b) = freqs(f)*2^(bands(b)/2) ;
        end
    end
    x = top_freqs(:) ;
    y = resps(:) ;
    f = d.top_freq_fit(neuron,:) ;
    x_hat = 3:0.1:20 ;
    y_hat =  f(1)*exp(-((sort(log(x_hat))-f(2))/f(3)).^2) ;
    ['top R_square = ' num2str(d.top_freq_R(neuron))]
    hold(Ax,'on')
    plot(Ax , x , y ,'.k')
    plot(Ax , x_hat , y_hat ,'r')
    set(Ax,'XScale','log','XTick',2.^(1:5))
    ylabel(Ax,'Firing rate (Hz)','FontSize',14)
    xlabel(Ax,'Frequency (kHz)','FontSize',14)
    xlim(Ax,[3 20])
    
end

function plot_bottom_freq_tuning(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder{handles.mouse_menu.Value} ;
    neuron = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;
    freqs = d.freqs ;
    bands = d.bandwidths ;
    resps = squeeze(d.amps_bands(neuron,:,2:end)) ;
    bottom_freqs = zeros(length(freqs),length(bands)) ;
    for f = 1 : length(freqs)
        for b = 1 : length(bands)
            bottom_freqs(f,b) = freqs(f)*2^(-bands(b)/2) ;
        end
    end
    x = bottom_freqs(:) ;
    y = resps(:) ;
    x_hat = 3:0.1:20 ;
    f = d.bottom_freq_fit(neuron,:) ;
    y_hat =  f(1)*exp(-((sort(log(x_hat))-f(2))/f(3)).^2) ;
    ['bottom R_square = ' num2str(d.bottom_freq_R(neuron))]
    hold(Ax,'on')
    plot(Ax , x , y ,'.k')
    plot(Ax , x_hat , y_hat ,'r')
    set(Ax,'XScale','log','XTick',2.^(1:5))
    ylabel(Ax,'Firing rate (Hz)','FontSize',14)
    xlabel(Ax,'Frequency (kHz)','FontSize',14)
    xlim(Ax,[3 20])
    
end

function plot_FBRA(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder{handles.mouse_menu.Value} ;
    neuron = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;

    r = squeeze(d.amps_bands(neuron,:,:))' ;

    imagesc(Ax , r)
    
    ylabel(Ax,'Bandwidth (% oct.)','FontSize',14)
    xlabel(Ax,'Frequency (kHz)','FontSize',14)
    set(Ax,'XTick', 1:size(r,2),'XTickLabel',strsplit(num2str(d.freqs,2)),...
        'YTick', 1:size(r,2),'YTickLabel',strsplit(num2str(100*[0 d.bandwidths],2)))
%     xlim(Ax,[0.5 11.5])
    
end

function plot_nFBRA(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder{handles.mouse_menu.Value} ;
    neuron = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;

    r = squeeze(d.amps_notches(neuron,:,2:end))' ;

    imagesc(Ax , r)
    
    ylabel(Ax,'Bandwidth (% oct.)','FontSize',14)
    xlabel(Ax,'Frequency (kHz)','FontSize',14)
    set(Ax,'XTick', 1:size(r,2),'XTickLabel',strsplit(num2str(d.freqs,2)),...
        'YTick', 1:size(r,2),'YTickLabel',strsplit(num2str(100*d.bandwidths,2)))
%     xlim(Ax,[0.5 11.5])
    
end

function plot_tuning_curve(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder{handles.mouse_menu.Value} ;
    neuron = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;
    
    x_hat = d.fit_parameters_for_disp(neuron,:) ;
    
    t = d.make_tuning_curve(x_hat(1),x_hat(2),x_hat(3),x_hat(4)) ;
    
    plot(Ax , d.sim_freqs , t)
    
%     axis(Ax,[4 16*2^0.5 -1.2 1.2])
    ylim(Ax,[-1.2 1.2])
    ylabel(Ax,'Weight (A.U.)','FontSize',14)
    xlabel(Ax,'Frequency (kHz)','FontSize',14)
    set(Ax , 'XScale' , 'log','XTick',2.^(2:6))
    
end

function plot_fitted_FBRA(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder{handles.mouse_menu.Value} ;
    neuron = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;
    
    x_hat = d.fit_parameters_for_disp(neuron,:) ;    
    r_hat = d.calc_sim_band_responses(x_hat(1),x_hat(2),x_hat(3),x_hat(4))' ;
        
    imagesc(Ax , r_hat)
    
    ylabel(Ax,'Bandwidth (% oct.)','FontSize',14)
    xlabel(Ax,'Frequency (kHz)','FontSize',14)
    set(Ax,'XTick', 1:size(r_hat,2),'XTickLabel',strsplit(num2str(d.freqs,2)),...
        'YTick', 1:size(r_hat,2),'YTickLabel',strsplit(num2str(100*[0 d.bandwidths],2)))
%     xlim(Ax,[0.5 11.5])
    
end

function plot_predicted_nFBRA(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder{handles.mouse_menu.Value} ;
    neuron = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;
    
    x_hat = d.fit_parameters_for_disp(neuron,:) ;    
    r_hat = d.calc_sim_notch_responses(x_hat(1),x_hat(2),x_hat(3),x_hat(4))' ;
        
    imagesc(Ax , r_hat)
    
    ylabel(Ax,'Bandwidth (% oct.)','FontSize',14)
    xlabel(Ax,'Frequency (kHz)','FontSize',14)
    set(Ax,'XTick', 1:size(r_hat,2),'XTickLabel',strsplit(num2str(d.freqs,2)),...
        'YTick', 1:size(r_hat,2),'YTickLabel',strsplit(num2str(100*d.bandwidths,2)))
%     xlim(Ax,[0.5 11.5])
    
end


function plot_fit_corrs_dist(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder ;
    
    hold(Ax,'on')
    
    corrs = [] ;
    shuff_corrs = [] ;

    for mouse = 1 : numel(d)
        dd = d{mouse} ;
        shuff_corrs = [shuff_corrs dd.fit_corrs_shuff(logical(dd.sig_cells_bands(1,:)))] ;
        corrs = [corrs dd.fit_corrs(logical(dd.sig_cells_bands(1,:)))] ;
    end
    
    histogram(Ax,corrs,0:0.1:1)
    histogram(Ax,shuff_corrs,0:0.1:1)
    ylim([0 90])
    set(Ax,'XTick',0:0.2:1,'YTick',[0 45 90])
    
    xlabel(Ax,'Corr. values','FontSize',14)
    ylabel(Ax,'# cells','FontSize',14)
    
    mean(corrs)
    median(corrs)
    std(corrs)
    length(corrs)

    mean(shuff_corrs)
    median(shuff_corrs)
    std(shuff_corrs)
    
end

function plot_notch_corrs_dist(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder ;
    
    hold(Ax,'on')
    
    corrs = [] ;
    shuff_corrs = [] ;

    for mouse = 1 : numel(d)
        dd = d{mouse} ;
        corrs = [corrs dd.notch_corrs(logical(dd.sig_cells_bands(1,:)))] ;
        shuff_corrs = [shuff_corrs dd.notch_corrs_shuff(logical(dd.sig_cells_bands(1,:)))] ;
    end
    

    histogram(Ax,corrs,-0.5:0.1:1)
    histogram(Ax,shuff_corrs,-0.5:0.1:1)
%         ylim([0 25])
    
    xlabel(Ax,'Corr. values','FontSize',14)
    ylabel(Ax,'# cells','FontSize',14)
    
    [~,p] = ttest(corrs,shuff_corrs)
    mean(corrs)
    median(corrs)
    std(corrs)
    length(corrs)
    
end

function plot_corrs_scatter(hObject,Ax)

    handles = guidata(hObject) ;
    d = handles.data_holder ;
    
    hold(Ax,'on')
    
    band_corrs = [] ;
    notch_corrs = [] ;
    
    handles.h_scat = {} ;
    for mouse = 1 : numel(d)
        dd = d{mouse} ;
        band_corrs = [band_corrs dd.fit_corrs(logical(dd.sig_cells_bands(1,:)))] ;
        notch_corrs = [notch_corrs dd.notch_corrs(logical(dd.sig_cells_bands(1,:)))] ;
        neurons = find(dd.sig_cells_bands(1,:)) ;
        for n = 1 : length(neurons)
            neur = neurons(n) ;
            handles.h_scat{end+1} = scatter(Ax,dd.fit_corrs(neur),dd.notch_corrs(neur),30,...
                'MarkerEdgeColor','k','MarkerFaceColor','k') ;
        end
    end
    
    xlabel(Ax,'Band corrs','FontSize',14)
    ylabel(Ax,'Notch corrs','FontSize',14)
    
    [temp,p] = corrcoef(band_corrs,notch_corrs) ;
    ['Corr. = ' num2str(temp(1,2))]
    ['p = ' num2str(p(1,2))]
    
    guidata(hObject,handles)
    
end

function plot_param_dist(hObject,Ax,i_param)

    handles = guidata(hObject) ;
    d = handles.data_holder ;
        
    param_vals = [] ;

    for mouse = 1 : numel(d)
        dd = d{mouse} ;
        x_hat = dd.fit_parameters(logical(dd.sig_cells_bands(1,:)) & logical(dd.fit_corrs > handles.min_corr),:) ;
        param_vals = [param_vals x_hat(:,i_param)'] ;
    end
    
    switch i_param
        case 1
            histogram(Ax,log(param_vals),log(2.^(2:0.4:6)))
            xlabel(Ax,'Mean values (kHz)','FontSize',14)
            set(Ax,'XTick',log(2.^(2:7)),'XTickLabel',strsplit(num2str(2.^(2:7))))
            ['Mean = ' num2str(mean(param_vals))]
            ['Standard deviation = ' num2str(std(param_vals))]
        case 2
            histogram(Ax,param_vals,0:0.1:1.2)
            xlabel(Ax,'Width values (oct.)','FontSize',14)
            ['Mean = ' num2str(mean(param_vals))]
            ['Standard deviation = ' num2str(std(param_vals))]
        case 3
            histogram(Ax,param_vals,-1:0.2:1)
            xlabel(Ax,'Skewness values','FontSize',14)
        case 4
            histogram(Ax,param_vals,0:0.1:1)
            xlabel(Ax,'I/E values','FontSize',14)
    end
    
    
end

function plot_skew_vs_IE(hObject,Ax)

    handles = guidata(hObject) ;
    hold(Ax,'on')
    
    d = handles.data_holder ;
        
    plot(Ax,[0 1],[0 0],'--k')
    plot(Ax,[0.5 0.5],[-1 1],'--k')
    
    handles.h_scat = {} ;
    for mouse = 1 : numel(d)
        dd = d{mouse} ;
        neurons = find(dd.sig_cells_bands(1,:)) ;
        for n = 1 : length(neurons)
            neur = neurons(n) ;
            x_hat = dd.fit_parameters(neur,:) ;
            if dd.fit_corrs(neur) > handles.min_corr
                handles.h_scat{end+1} = scatter(Ax,x_hat(4),x_hat(3),30,...
                    'MarkerEdgeColor','k','MarkerFaceColor','k') ;
            end
            
        end
    end
    
    set(Ax,'XTick',[0 0.5 1],'YTick',[-1 0 1])
    xlabel(Ax,'I/E','FontSize',14)
    ylabel(Ax,'Skewness','FontSize',14)
    
    guidata(hObject,handles)
    
end

function plot_skew_vs_IE_quad_dist(hObject,Ax)

    handles = guidata(hObject) ;
    hold(Ax,'on')
    
    d = handles.data_holder ;
        
    Skew = [] ;
    IE = [] ;
    
    for mouse = 1 : numel(d)
        dd = d{mouse} ;
        neurons = find(dd.sig_cells_bands(1,:)) ;
        for n = 1 : length(neurons)
            neur = neurons(n) ;
            if  dd.fit_corrs(neur) >= handles.min_corr
                x_hat = dd.fit_parameters(neur,:) ;
                Skew(end+1) = x_hat(3) ;
                IE(end+1) = x_hat(4) ;
            end
        end
    end
    
    q(1) = 100*sum(Skew > 0 & IE > 0.5)/length(Skew) ;
    q(2) = 100*sum(Skew > 0 & IE < 0.5)/length(Skew) ;
    q(3) = 100*sum(Skew < 0 & IE < 0.5)/length(Skew) ;
    q(4) = 100*sum(Skew < 0 & IE > 0.5)/length(Skew) ;
    
    bar(Ax,1:4,q)
    set(Ax,'XTick',1:4,'XTickLabel',{'Q1','Q2','Q3','Q4'},'FontSize',14)
    ylabel(Ax,'% of cells','FontSize',14)
    
    guidata(hObject,handles)
    
end

function plot_dR_fit_per_quad(hObject,Ax)

    handles = guidata(hObject) ;
    hold(Ax,'on')
    
    d = handles.data_holder ;
        
    top_vs_bottom = [] ;
    Skew = [] ;
    IE = [] ;
    
    for mouse = 1 : numel(d)
        dd = d{mouse} ;
        neurons = find(dd.sig_cells_bands(1,:)) ;
        for n = 1 : length(neurons)
            neur = neurons(n) ;
            if  dd.fit_corrs(neur) >= handles.min_corr
                x_hat = dd.fit_parameters(neur,:) ;
                top_vs_bottom(end+1) = dd.top_freq_R(neur)-dd.bottom_freq_R(neur);
                Skew(end+1) = x_hat(3) ;
                IE(end+1) = x_hat(4) ;
            end
        end
    end
    
    q1 = top_vs_bottom(Skew > 0 & IE > 0.5) ;
    q_else = top_vs_bottom(~(Skew > 0 & IE > 0.5)) ;
    
    bar(Ax,1:2,[mean(q1) mean(q_else)])
    errorbar(Ax,1:2,[mean(q1) mean(q_else)],[std(q1)/sqrt(length(q1)) std(q_else)/sqrt(length(q_else))],'.k')
    set(Ax,'YTick',[0 0.3],'XTick',1:2,'XTickLabel',{'Q1','Q2-4'},'FontSize',14)
    ylabel(Ax,'top-bottom fit','FontSize',14)
    ylim(Ax,[-0.1 0.3])
    
    p = ranksum(q1,q_else)
    
    guidata(hObject,handles)
    
end

function plot_top_freq_vs_bottom_freq(hObject,Ax)

    handles = guidata(hObject) ;
    hold(Ax,'on')
    
    d = handles.data_holder ;
        
    top_freq_R = [] ;
    bottom_freq_R = [] ;
    
    handles.h_scat = {} ;
    for mouse = 1 : numel(d)
        dd = d{mouse} ;
        neurons = find(dd.sig_cells_bands(1,:)) ;
        for n = 1 : length(neurons)
            neur = neurons(n) ;
            temp = dd.fit_parameters(neur,1) ;
            if  dd.fit_corrs(neur) >= handles.min_corr
                x_hat = dd.fit_parameters(neur,:) ;
                top_freq_R(end+1) = dd.top_freq_R(neur) ;
                bottom_freq_R(end+1) = dd.bottom_freq_R(neur) ;
                if x_hat(3) > 0 && x_hat(4) > 0.5
                    handles.h_scat{end+1} = scatter(Ax,bottom_freq_R(end),top_freq_R(end),30,...
                        'MarkerEdgeColor','k','MarkerFaceColor','k') ;
                else
                    handles.h_scat{end+1} = scatter(Ax,bottom_freq_R(end),top_freq_R(end),30,...
                        'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1]) ;
                end
                set(handles.h_scat{end},'ButtonDownFcn',{@scatter_func , hObject , mouse , n}) ;
            end
        end
    end
    
    
    LIMS = [0 1] ;
    plot(Ax,LIMS,LIMS,'r')
    ylabel('Top freq. R','FontSize',14)
    xlabel('Bottom freq. R','FontSize',14)
    set(Ax,'XTick',[0 0.5 1],'YTick',[0 0.5 1])
   
    guidata(hObject,handles)
    
end

% ========================== Helper functions ============================

function band = create_band(freq,Bandwidth)

    % inputs:
    % freq (middle frequency kHz)
    % Bandwidth (fraction of an octave)

    
    unit_amp = 0.05 ;   % volts
    band_res = 0.001 ;  % fraction oct.
    Fs = 500 ;          % kHz
    dur = 100 ;         % ms

    band = zeros(1,Fs*dur) ;
    for i = (-Bandwidth/2):band_res:(Bandwidth/2)
        band = band + unit_amp * sin( (rand() + freq*2^i / Fs * (1:dur*Fs)) * 2*pi ) ;
    end
    
end

function scroll_func(hObject,eventdata)

    handles = guidata(hObject) ;
    
    fig_pos = get(hObject,'Position') ;
    mouse_pos = get(hObject,'CurrentPoint') ;
    
    anlss_type_menu_pos(1) = fig_pos(3) * handles.anlss_type_menu.Position(1) ;
    anlss_type_menu_pos(2) = fig_pos(3) * (handles.anlss_type_menu.Position(1) + handles.anlss_type_menu.Position(3)) ;
    anlss_type_menu_pos(3) = fig_pos(4) * handles.anlss_type_menu.Position(2) ;
    anlss_type_menu_pos(4) = fig_pos(4) * (handles.anlss_type_menu.Position(2) + handles.anlss_type_menu.Position(4)) ;
    
    neuron_menu_pos(1) = fig_pos(3) * handles.neuron_menu.Position(1) ;
    neuron_menu_pos(2) = fig_pos(3) * (handles.neuron_menu.Position(1) + handles.neuron_menu.Position(3)) ;
    neuron_menu_pos(3) = fig_pos(4) * handles.neuron_menu.Position(2) ;
    neuron_menu_pos(4) = fig_pos(4) * (handles.neuron_menu.Position(2) + handles.neuron_menu.Position(4)) ;
    
    mouse_menu_pos(1) = fig_pos(3) * handles.mouse_menu.Position(1) ;
    mouse_menu_pos(2) = fig_pos(3) * (handles.mouse_menu.Position(1) + handles.mouse_menu.Position(3)) ;
    mouse_menu_pos(3) = fig_pos(4) * handles.mouse_menu.Position(2) ;
    mouse_menu_pos(4) = fig_pos(4) * (handles.mouse_menu.Position(2) + handles.mouse_menu.Position(4)) ;
    
    if (mouse_pos(1) > anlss_type_menu_pos(1)) && (mouse_pos(1) < anlss_type_menu_pos(2)) &&...
       (mouse_pos(2) > anlss_type_menu_pos(3)) && (mouse_pos(2) < anlss_type_menu_pos(4))
   
        val = handles.anlss_type_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.anlss_type_menu.String) ;
   
        handles.anlss_type_menu.Value = min([max([new_val 1]) max_val]) ;
        
        anlss_type_menu_Callback(hObject,0,0)
        
    elseif (mouse_pos(1) > neuron_menu_pos(1)) && (mouse_pos(1) < neuron_menu_pos(2)) &&...
           (mouse_pos(2) > neuron_menu_pos(3)) && (mouse_pos(2) < neuron_menu_pos(4))
   
        val = handles.neuron_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.neuron_menu.String) ;
   
        handles.neuron_menu.Value = min([max([new_val 1]) max_val]) ;
        
        neuron_menu_Callback(hObject,0,0)
        
    elseif (mouse_pos(1) > mouse_menu_pos(1)) && (mouse_pos(1) < mouse_menu_pos(2)) &&...
           (mouse_pos(2) > mouse_menu_pos(3)) && (mouse_pos(2) < mouse_menu_pos(4))
   
        val = handles.mouse_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.mouse_menu.String) ;
   
        handles.mouse_menu.Value = min([max([new_val 1]) max_val]) ;
        
        mouse_menu_Callback(hObject,0,0)
                       
    end
    
end

function disp_parameters(hObject,x_hat,fit_corr,notch_corr)

    handles = guidata(hObject) ;
    d = handles.data_holder{handles.mouse_menu.Value} ;
    
    handles.corr_txt.String(8:end) = [] ;
    handles.corr_txt.String = [handles.corr_txt.String num2str(fit_corr,2)] ;
    
    handles.mean_txt.String(8:end) = [] ;
    handles.mean_txt.String = [handles.mean_txt.String num2str(x_hat(1),2) ' kHz'] ;
    
    handles.width_txt.String(9:end) = [] ;
    handles.width_txt.String = [handles.width_txt.String num2str(x_hat(2),2) ' oct.'] ;
    
    handles.skew_txt.String(8:end) = [] ;
    handles.skew_txt.String = [handles.skew_txt.String num2str(x_hat(3),2)] ;
    
    handles.IE_txt.String(7:end) = [] ;
    handles.IE_txt.String = [handles.IE_txt.String num2str(x_hat(4),2)] ;
    
    handles.notch_corr_txt.String(14:end) = [] ;
    handles.notch_corr_txt.String = [handles.notch_corr_txt.String num2str(notch_corr,2)] ;
    
end


% ========================== Callbacks ===================================

function load_data_button_Callback(hObject, eventdata, handles)

    [file,path] = uigetfile ;
    
    if file
        handles.load_data_button.Enable = 'off' ;
        load([path file],'d')
        handles.data_holder = d ;
        
        handles.mouse_menu.String = {} ;
        for i = 1 : numel(d)
            handles.mouse_menu.String{end+1} = d{i}.mouse_name ;
        end
        handles.mouse_menu.Value = 1 ;
        guidata(hObject,handles)
        initialize_all(hObject)
        guidata(hObject,handles)
    end
    
end

function mouse_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    
    d = handles.data_holder{handles.mouse_menu.Value} ;
    
    handles.neuron_menu.String = {} ;
    for i = 1 : size(d.sig_cells_bands,2)    
        if d.sig_cells_bands(1,i)     
            handles.neuron_menu.String{end+1} = ['Neuron ' num2str(i)] ;
        end       
    end
    if isempty(handles.neuron_menu.String)
        handles.neuron_menu.String{1} = 'No neurons' ;
    end
    
    handles.neuron_menu.Value = 1 ;    
    neuron_menu_Callback(hObject, 0 , 0)
        
end

function neuron_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;

    d = handles.data_holder{handles.mouse_menu.Value} ;
    if ~strcmp(handles.neuron_menu.String,'No neurons')
        neuron = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;
    else
        return
    end
    
    temp = string(d.area_acronyms) ;
    handles.area_name_txt.String = temp{neuron} ;
    anlss_type_menu_Callback(hObject,0,0)
    
    disp_parameters(hObject , d.fit_parameters(neuron,:) , d.fit_corrs(neuron) , d.notch_corrs(neuron))
    
    clc


end

% ========================== Creation functions ==========================

function mouse_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function neuron_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

end

function anlss_type_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end
