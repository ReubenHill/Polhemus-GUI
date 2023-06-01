function choice = doubletapdialog(currentval)

    d = dialog('Position',[300 300 260 150],'Name','Double tap');

    slider = uicontrol('Parent',d,...
           'Style','slider',...
           'Position',[31 60 200 25],...
           'Min',0,...
           'Max',1,...
           'Value', currentval,...
           'Callback',@slider_callback);

    txt = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[31 100 210 30],...
           'String', sprintf('Use slider to set warning distance.\nCurrent value: %0.2g cm.', currentval));

    fun = @(~,e)set(txt,'String',sprintf('Use slider to set warning distance.\nCurrent value: %0.2g cm.', get(e.AffectedObject,'Value')));

    addlistener(slider, 'Value', 'PostSet', fun);

    btn = uicontrol('Parent',d,...
           'Position',[99 20 70 25],...
           'String','Close',...
           'Callback','delete(gcf)');

    choice = 'Red';

    % Wait for d to close before running to completion
    uiwait(d);

    function slider_callback(slider,event)
        choice = slider.Value;
    end

end