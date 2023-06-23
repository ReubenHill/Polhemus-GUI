function [choice, enabled] = doubletapdialog(currentchoice, currentenabled)

    d = dialog('Position',[300 300 260 180],'Name','Double tap');

    checkbox = uicontrol('Parent', d, ...
            'Style', 'checkbox', ...
            'String', 'Enable double tap error?', ...
            'Position', [68 70 200 170], ...
            'Value', currentenabled, ...
            'Callback', @checkbox_callback);

    slider = uicontrol('Parent',d,...
           'Style','slider',...
           'Position',[31 60 200 25],...
           'Min',0,...
           'Max',1,...
           'Value', currentchoice,...
           'Callback',@slider_callback);

    txt = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[31 100 210 30],...
           'String', sprintf('Use slider to set error distance.\nCurrent value: %0.2g cm.', currentchoice));

    fun = @(~,e)set(txt,'String',sprintf('Use slider to set error distance.\nCurrent value: %0.2g cm.', get(e.AffectedObject,'Value')));

    addlistener(slider, 'Value', 'PostSet', fun);

    btn = uicontrol('Parent',d,...
           'Position',[99 20 70 25],...
           'String','Close',...
           'Callback','delete(gcf)');

    choice = currentchoice;
    enabled = currentenabled;

    if currentenabled
        set(txt,'Enable','on');
        set(slider,'Enable','on');
    else
        set(txt,'Enable','off');
        set(slider,'Enable','off');
    end

    % Wait for d to close before running to completion
    uiwait(d);


    function slider_callback(slider,event)
        choice = slider.Value;
    end

    function checkbox_callback(checkbox,event)
        enabled = checkbox.Value;
        if enabled
            set(txt,'Enable','on');
            set(slider,'Enable','on');
        else
            set(txt,'Enable','off');
            set(slider,'Enable','off');
        end
    end

end