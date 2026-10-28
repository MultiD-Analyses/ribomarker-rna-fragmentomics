/* Settings above the plots of the technical QC report (qcReport.py inlines this file into the report).
   Martin 2026-10-26: settings (Reads / Percentage, Linear / Log scale, % of natural reads / % within class)
   are shown the same way everywhere, on a line of their own above the plot's other buttons, and look different
   from them. MultiQC's own toggles (percentage, log scale) are hidden and driven by these buttons; plots with a
   mode setting (a div.rm-modes in their description) have their datasets as mode x choice, and MultiQC's dataset
   buttons are replaced by one row of choices. Tables in a div.rm-after-plot are moved below their plot. */
(function () {
  function settingsLine(container) {
    var line = container.querySelector(':scope > .rm-settings-line');
    if (!line) {
      line = document.createElement('div');
      line.className = 'rm-settings-line';
      container.insertBefore(line, container.firstChild);
    }
    return line;
  }

  function buttonRow(cls, options, active, onPick) {
    var group = document.createElement('div');
    group.className = 'btn-group ' + cls;
    var buttons = options.map(function (text, i) {
      var b = document.createElement('button');
      b.type = 'button';
      b.textContent = text;
      b.className = cls === 'rm-setting-group' ? 'btn btn-sm rm-setting' : 'btn btn-outline-secondary btn-sm';
      b.addEventListener('click', function (e) {
        e.preventDefault();
        onPick(i);
        buttons.forEach(function (x, j) { x.classList.toggle('active', j === i); });
      });
      group.appendChild(b);
      return b;
    });
    buttons[active].classList.add('active');
    return group;
  }

  function setting(container, label, options, active, onPick) {
    var wrap = document.createElement('span');
    wrap.className = 'rm-setting-wrap';
    var name = document.createElement('span');
    name.className = 'rm-setting-label';
    name.textContent = label + ':';
    wrap.appendChild(name);
    wrap.appendChild(buttonRow('rm-setting-group', options, active, onPick));
    settingsLine(container).appendChild(wrap);
  }

  document.addEventListener('DOMContentLoaded', function () {
    // MultiQC's toggles: percentage and log scale (RiboMarker plots and the tool sections alike)
    document.querySelectorAll('button.interactive-switch-group.percent-switch, button.interactive-switch-group.log10-switch')
      .forEach(function (toggle) {
        var percent = toggle.classList.contains('percent-switch');
        var ours = String(toggle.dataset.plotAnchor || '').indexOf('ribomarker_') === 0;
        var options = percent ? [ours ? 'Reads' : 'Counts', 'Percentage'] : ['Linear', 'Log'];
        toggle.style.display = 'none';
        setting(toggle.closest('.col-12') || toggle.parentNode, percent ? 'Values' : 'Scale', options,
                toggle.classList.contains('active') ? 1 : 0, function (i) {
                  if ((i === 1) !== toggle.classList.contains('active')) { toggle.click(); }
                });
      });
    // Tables collapsed below a plot: MultiQC puts a section's description above its plot, so they are moved
    // below it (div.rm-after-plot)
    document.querySelectorAll('div.rm-after-plot').forEach(function (div) {
      var section = div.closest('.mqc-section');
      if (section) { section.appendChild(div); }
    });
    // A mode setting over MultiQC datasets (mode x choice)
    document.querySelectorAll('div.rm-modes').forEach(function (cfg) {
      var anchor = cfg.dataset.anchor;
      var modes = JSON.parse(cfg.dataset.modes);
      var native = Array.prototype.slice.call(document.querySelectorAll(
        '.dataset-switch-group button[data-plot-anchor="' + anchor + '"]'));
      if (!native.length) { return; }
      if (native.length % modes.length) { return; }   // not mode x choice: leave MultiQC's buttons as they are
      var n = native.length / modes.length;
      var group = native[0].parentNode;
      var container = group.closest('.col-12') || group.parentNode;
      var state = {mode: 0, choice: 0};
      native.forEach(function (b, i) { if (b.classList.contains('active')) { state = {mode: Math.floor(i / n), choice: i % n}; } });
      function show() {
        var target = native[state.mode * n + state.choice];
        if (!target.classList.contains('active')) { target.click(); }
      }
      var choices = native.slice(0, n).map(function (b) { return b.textContent; });
      group.style.display = 'none';
      container.insertBefore(buttonRow('rm-choice-group', choices, state.choice, function (i) { state.choice = i; show(); }),
                             group);
      setting(container, cfg.dataset.label || 'Values', modes, state.mode, function (i) { state.mode = i; show(); });
    });
  });

  // Colour samples by a sample-sheet column (Martin 2026-10-26): the chosen column's values become MultiQC
  // highlight patterns, which colour the samples in every plot and table; plots whose series are the samples
  // (series_plots) get the colours directly. Waits until MultiQC has made its plots and toolbox.
  function escapeRegex(text) { return text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); }

  // MultiQC's bar plots take a category colour as "r,g,b" (bar.js wraps it in rgb()); a hex colour there is invalid
  function rgbList(hex) {
    var n = parseInt(hex.slice(1), 16);
    return [(n >> 16) & 255, (n >> 8) & 255, n & 255].join(',');
  }

  var chosen = {};                              // column -> {value: colour picked in the key}; kept while the page is open

  function colourSamples(info, column, key) {
    var colourOf = {};                          // sample -> colour; {} = no colours
    var groups = [];
    key.textContent = '';
    if (column) {
      var values = info.columns[column];
      groups = values.filter(function (v, i) { return values.indexOf(v) === i; });
      if (groups.length > info.palette.length) {
        key.textContent = groups.length + ' different values: too many to colour (at most ' + info.palette.length + ').';
        groups = [];
      } else {
        chosen[column] = chosen[column] || {};
        var colourOfGroup = {};
        groups.forEach(function (g, i) { colourOfGroup[g] = chosen[column][g] || info.palette[i]; });
        info.samples.forEach(function (s, i) { colourOf[s] = colourOfGroup[values[i]]; });
        groups.forEach(function (g) {           // the key: one colour picker per value
          var item = document.createElement('label');
          item.className = 'rm-colour-chip';
          item.style.color = colourOfGroup[g];
          item.title = 'Click the square to choose the colour of ' + g;
          var picker = document.createElement('input');
          picker.type = 'color';
          picker.value = colourOfGroup[g];
          picker.addEventListener('change', function () {
            chosen[column][g] = picker.value;
            colourSamples(info, column, key);
          });
          item.appendChild(picker);
          item.appendChild(document.createTextNode(' ' + g));
          key.appendChild(item);
        });
      }
    }
    var filters = $('#mqc_col_filters').empty();
    groups.forEach(function (g) {
      var names = info.samples.filter(function (s, i) { return info.columns[column][i] === g; }).map(escapeRegex);
      filters.append(window.make_colorsamples_filter('^(' + names.join('|') + ')$', colourOf[info.samples[info.columns[column].indexOf(g)]]));
    });
    $('#mqc_cols .mqc_regex_mode input').prop('checked', true);
    info.series_plots.forEach(function (anchor) {
      var plot = window.mqc_plots[anchor];
      if (!plot) { return; }
      plot.datasets.forEach(function (dataset) {
        (dataset.cats || []).forEach(function (cat) {
          if (cat.rmOriginal === undefined) { cat.rmOriginal = cat.color; }
          cat.color = colourOf[cat.name] ? rgbList(colourOf[cat.name]) : cat.rmOriginal;
        });
      });
    });
    window.apply_mqc_highlights();               // re-renders every plot and table
  }

  function initColourSelector(tries) {
    var box = document.querySelector('div.rm-colour-by');
    if (!box) { return; }
    if (!(window.mqc_plots && window.apply_mqc_highlights && window.make_colorsamples_filter)) {
      if (tries < 200) {
        setTimeout(function () { initColourSelector(tries + 1); }, 100);
      } else {                                  // MultiQC's report scripts changed (written for MultiQC 1.35)
        box.textContent = 'Colour samples by: not available in this report (it needs the scripts of MultiQC 1.35); '
          + 'the toolbox\u2019s Highlight tab can colour samples instead.';
      }
      return;
    }
    var bytes = Uint8Array.from(atob(box.dataset.config), function (c) { return c.charCodeAt(0); });
    var info = JSON.parse(new TextDecoder('utf-8').decode(bytes));
    var label = document.createElement('label');
    label.className = 'rm-setting-label';
    label.textContent = 'Colour samples by:';
    var select = document.createElement('select');
    select.className = 'form-select form-select-sm rm-colour-select';
    ['(none)'].concat(Object.keys(info.columns)).forEach(function (name, i) {
      var option = document.createElement('option');
      option.value = i === 0 ? '' : name;
      option.textContent = name;
      select.appendChild(option);
    });
    select.value = info.first;
    var key = document.createElement('span');
    key.className = 'rm-colour-key';
    box.appendChild(label);
    box.appendChild(select);
    box.appendChild(key);
    select.addEventListener('change', function () { colourSamples(info, select.value, key); });
    colourSamples(info, select.value, key);
  }
  // after MultiQC's own start-up (jQuery ready handlers: tables, toolbox, saved settings), so that the tables
  // listen already when the first colours are applied
  document.addEventListener('DOMContentLoaded', function () {
    $(function () { setTimeout(function () { initColourSelector(0); }, 0); });
  });
})();
