var app = app || {};

(function ($) {
	/***************
	 * OverView.js *
	 ***************/
	app.OverView = app.BaseView.extend({
		// Our template for the line of statistics at the bottom of the app.
		template: _.template($('#overview-template').html()),

		// Delegated events for creating new items, and clearing completed ones.
		events: {
			'click #mesmer-item a': 'headClick',
			'dblclick #mesmer-title': 'rename',
			'blur #mesmer-title input': 'updateTitle',
			'keydown #mesmer-title input': 'onEnter'
		},

		render: function () {
            this.$el.appendTo($('body'));
            this.$el.addClass('container');
			this.$el.html(this.template(this.model.toJSON()));
			this.$('#mesmer-title .list-group-item').tooltip({
				trigger: 'hover focus',
				title: 'double click to edit title',
				placement: 'top'
			});
			return this;
		},

		initialize: function () {
			this.moleculeList = this.model.get('moleculeList');
			this.reactionList = this.model.get('reactionList');
			this.control = this.model.get('control');
			this.conditions = this.model.get('conditions');
			this.modelParameters = this.model.get('modelParameters');
			
			this.listenTo(this.model, 'change', this.myRender);
			this.listenTo(this.model, 'destroy', this.remove);

			this.reactionListView = new app.ReactionListView(this.reactionList, this.moleculeList);
			this.controlView = new app.ControlView({model: this.control});
			this.conditionsView = new app.ConditionsView({model: this.conditions}, this.moleculeList);
			this.moleculeListView = new app.MoleculeListView(this.moleculeList);
			this.modelParametersView = new app.ModelParametersView({model: this.modelParameters});

			this.render().$('#content').
			append(this.reactionListView.render().el).
			append(this.controlView.render().el).
			append(this.conditionsView.render().el).
			append(this.modelParametersView.render().el).
			append(this.moleculeListView.render().el);

			this.reactionListView.$el.hide();
            this.controlView.$el.hide();
            this.conditionsView.$el.hide();
            this.modelParametersView.$el.hide();
            this.$('#molecule-list').addClass('active');
		},

		rename: function (e) {
            var el = $(e.currentTarget);
			el.find('.list-group-item').hide();
			el.find('input').show().focus();
		},

		headClick: function(e){
			var el = $(e.currentTarget);
			this.$('#mesmer-item .active').removeClass('active');
			el.addClass('active');

			this.$('#content').children().hide();
			switch (el.attr('id')){
				case 'molecule-list':
					this.moleculeListView.render();
					this.moleculeListView.show();
					break;

				case 'reaction-list':
					this.reactionListView.render();
					this.reactionListView.show();
					break;

				case 'conditions':
					this.conditionsView.render();
					this.conditionsView.show();
					break;

				case 'model-parameters':
					this.modelParametersView.render();
					this.modelParametersView.show();
					break;

				case 'control':
					this.controlView.render();
					this.controlView.show();
					break;
			}
		},

		updateTitle: function(e){
            var $e = $(e.currentTarget);
            if (!this.checkFilename($e.val())){
                $e.popover({
                    content: "Title cannot contains characters like <, >, /, \\, |, :, \", *, ?, !",
                    placement: "bottom",
                    trigger: "manual"
                });
                $e.popover("show");
                setTimeout(function(){$e.popover("hide")}, 2000);
            } else {
                this.inputDone(e);
			    $e.hide().prev().show().text(this.model.get('title'));
            }
		},

		onEnter: function(e){
			if (e.which === 13){
				this.updateTitle(e);
			}
		},

        myRender: function(){
            this.$('#mesmer-title .list-group-item').text(this.model.get('title'));
            this.$('#mesmer-title input').val(this.model.get('title'))
            $('title').text(this.model.get('title'));
        },

        checkFilename: function(filename) {
            var forbidden = "<>/\\|:\"*?";
            if (filename.match(/<|>|\/|\\|:|\||"|\*|\?/)){
                return false;
            } else {
                return true;
            }
        }
	});

	/***********************
	 * MoleculeListView.js *
	 ***********************/
	app.MoleculeListView = app.BaseView.extend({

		template: _.template($('#molecule-list-property-template').html()),

		events: {
			'click button[action ^= "new-molecule"]': 'newMolecule',
            'click ul[target="library"] li a': 'addFromLibrary'
		},

		initialize: function(moleculeList) {
			this.moleculeList = moleculeList;
			this.moleculeView = new app.MoleculeView({});
			this.listenTo(this.moleculeList, 'add', this.addMolecule);
			this.listenTo(this.moleculeList, 'select', this.selectMolecule);
			this.listenTo(this.moleculeList, 'remove', this.selectAnother);
            this.listenTo(this.moleculeList, 'reset', this.reset);
            this.listenTo(this.moleculeList, 'library', this.showLibrary);
		},

		render: function() {
			this.$el.html(this.template());
			
			var molecules = this.moleculeList.models;
			for (i = 0; i < molecules.length; i++){
				var view = new app.MoleculeItemView({model: molecules[i]});
				this.$('#molecule-list-area').append(view.render().el);
			}
			if (i === 0){
				this.noMolecule();
			}
			else{
				this.moleculeList.select(this.moleculeList.selected || molecules[0]);
			}

            if (this.moleculeList.library)
                this.showLibrary();
            return this;
		},

        showLibrary: function(){
            var ul = this.$('ul[target^="library"]');
            for (key in this.moleculeList.library){
                var li = '<li><a href="#" mol-id="' + key + '">' + key + '</a></li>';
                ul.append(li);
            }
        },

        reset: function(){
            if (this.moleculeList.length > 0){
                this.hasMolecule();
                this.render();
            } else {
                this.noMolecule();
            }
        },

		newMolecule: function(e) {
			var molecule = new app.Molecule({id: 'M' + this.moleculeList.seq++});
			this.moleculeList.add(molecule);
			this.moleculeList.select(molecule);
			e.stopPropagation();
		},

        addFromLibrary: function(e) {
            var molID = $(e.currentTarget).attr('mol-id');
            var molecule = new app.Molecule({id: molID});
            var molFromLibrary = this.moleculeList.library[molID];
            for (key in molFromLibrary){
                if (molecule.get(key) != undefined)
                    molecule.set(key, molFromLibrary[key].value || molFromLibrary[key]);
                if (key == "ZPE")
                    molecule.set('ZPE_unit', molFromLibrary[key].units || "cm-1");
            }
            var mol = this.moleculeList.add(molecule);
			this.moleculeList.select(mol);
        },

		addMolecule: function(molecule) {
			var itemView = new app.MoleculeItemView({model: molecule});
			this.moleculeList.models.length === 1 && this.hasMolecule();
			this.$('#molecule-list-area').append(itemView.render().el);
		},

		selectMolecule: function(molecule){
			this.moleculeView.setModel(molecule);
			this.$('#molecule-list-area button').removeClass('active');
			this.$('#molecule-property').html(this.moleculeView.el);
			this.moleculeView.delegateEvents();
		},

		selectAnother: function(molecule){
			if (this.moleculeList.length === 0){
				this.noMolecule();
			}
			else {
				this.moleculeList.selected === molecule && 
				this.moleculeList.select(this.moleculeList.models[0]);
			}
		},

		noMolecule: function(){
			this.$('#no-molecule-message').slideDown('fast');
			this.$('#molecule-property').hide();
			this.$('#molecule-list-area').hide();
		},

		hasMolecule: function(){
			this.$('#molecule-property').show();
			this.$('#molecule-list-area').show();
			this.$('#no-molecule-message').hide();
		},

		show: function() {
			this.$el.show();
		}
	});

	/***********************
	 * MoleculeItemView.js *
	 ***********************/
	app.MoleculeItemView = app.BaseView.extend({
		tagName: 'span',
		
		template: _.template($('#molecule-item-template').html()),

		events: {
			'click button[action ^= "select"]': 'select',
			'click button[action ^= "delete"]': 'del'
		},

		initialize: function() {
			this.model.view = this;
		
			this.listenTo(this.model, 'change', this.render);
			this.listenTo(this.model, 'destroy', this.remove);
			this.listenTo(this.model, 'selected', this.setSelected);
			
		},

		render: function() {
			this.$el.html(this.template(this.model.toJSON()));
			this.model === this.model.collection.selected && this.setSelected();
			return this;
		},

		select: function(e) {
			this.model.collection.select(this.model);
		},

		del: function(e) {
			this.model.destroy();
		},
		
		setSelected: function() {
			this.$('button[action ^= "select"]').addClass("active");
		}
	});

	/*******************
	 * MoleculeView.js *
	 *******************/
	app.MoleculeView = app.BaseView.extend({
		template: _.template($('#molecule-property-template').html()),

		events: {
			'click ul[role ^= "menu"] a': 'selectMenu',
			'click button[data-button ^= "data-button"]': 'buttonClick',
			'blur input': 'inputDone',
			'keydown input': 'onEnter'
		},


		initialize: function(o, moleculeList) {
			this.moleculeList = moleculeList;
			this.model && this.listenTo(this.model, 'change', this.render);
		},


		render: function() {
			this.$el.html(this.template(this.model.toJSON()));
			this.hideSomething();
			return this;
		},

		hideSomething: function() {
			var items = app.constant.moleculeParameter[this.model.get('type')];
			for (i = 0; i < items.length; i++){
				this.$("[target = " + items[i] + "]").parent().hide();
			}
            if (this.model.get('type') == 'sink' || this.model.get('type') == 'bathGas'){
                this.$("[target=DOSCMethod]").parent().parent().hide();
            }
		},

		setModel: function (model) {
			this.stopListening();
			this.model = model;
			this.listenTo(this.model, 'change', this.render);
			this.render();
		}
	});

	/***********************
	 * ReactionListView.js *
	 ***********************/
	app.ReactionListView = app.BaseView.extend({

		template: _.template($('#reaction-list-property-template').html()),

		events: {
			'click button[action ^= "new-reaction"]': 'newReaction'
		},

		initialize: function(reactionList, moleculeList) {
			this.reactionList = reactionList;
			this.moleculeList = moleculeList;
			this.reactionView = new app.ReactionView({}, this.moleculeList);
			this.listenTo(this.reactionList, 'add', this.addReaction);
			this.listenTo(this.reactionList, 'select', this.selectReaction);
			this.listenTo(this.reactionList, 'remove', this.selectAnother);
            this.listenTo(this.reactionList, 'reset', this.reset);
		},

		render: function() {
			this.$el.html(this.template());
			
			var reactions = this.reactionList.models;
			for (i = 0; i < reactions.length; i++){
				var view = new app.ReactionItemView({model: reactions[i]});
				this.$('#reaction-list-area').append(view.render().el);
			}
			if (i === 0){
				this.noReaction();
			}
			else{
				this.reactionList.select(this.reactionList.selected || reactions[0]);
			}
			
			return this;
		},

        reset: function(){
            if (this.reactionList.length > 0){
                this.hasReaction();
                this.render();
            } else {
                this.noReaction();
            }
        },

		newReaction: function(e) {
			var reaction = new app.Reaction({id: 'R' + this.reactionList.seq++});
			this.reactionList.add(reaction);
			this.reactionList.select(reaction);
			e.stopPropagation();
		},

		addReaction: function(reaction) {
			var itemView = new app.ReactionItemView({model: reaction});
			this.reactionList.models.length === 1 && this.hasReaction();
			this.$('#reaction-list-area').append(itemView.render().el);
		},

		selectReaction: function(reaction){
			this.reactionView.setModel(reaction);
			this.$('#reaction-list-area button').removeClass('active');
			this.$('#reaction-property').html(this.reactionView.el);
			this.reactionView.delegateEvents();
		},

		selectAnother: function(reaction){
			if (this.reactionList.length === 0){
				this.noReaction();
			}
			else {
				this.reactionList.selected === reaction && 
				this.reactionList.select(this.reactionList.models[0]);
			}
		},

		noReaction: function(){
			this.$('#no-reaction-message').slideDown('fast');
			this.$('#reaction-property').hide();
			this.$('#reaction-list-area').hide();
		},

		hasReaction: function(){
			this.$('#reaction-property').show();
			this.$('#reaction-list-area').show();
			this.$('#no-reaction-message').hide();
		},

		show: function() {
			this.$el.show();
		}
	});
	/*************************
	 *  ReactionItemView.js  *
	 *************************/
	 
	app.ReactionItemView = app.BaseView.extend({
		tagName: 'span',
		
		template: _.template($('#reaction-item-template').html()),

		events: {
			'click button[action ^= "select"]': 'select',
			'click button[action ^= "delete"]': 'del'
		},

		initialize: function() {
			this.model.view = this;
			
			this.listenTo(this.model, 'change', this.render);
			this.listenTo(this.model, 'destroy', this.remove);
			this.listenTo(this.model, 'selected', this.setSelected);
			
			//this.$el.addClass('list-group-item list-group-item-sm');
		},

		// Re-render the titles of the todo item.
		render: function() {
			this.$el.html(this.template(this.model.toJSON()));
			this.model === this.model.collection.selected && this.setSelected();
			return this;
		},

		select: function(e) {
			this.model.collection.select(this.model);
		},

		del: function(e) {
			//this.model.collection.trigger('
			this.model.destroy();
		},
		
		setSelected: function() {
			this.$('button[action ^= "select"]').addClass("active");
		}
	});

	/*******************
	 * ReactionView.js *
	 *******************/
	app.ReactionView = app.BaseView.extend({
		template: _.template($('#reaction-property-template').html()),

		events: {
			'click ul[role ^= "menu"] a': 'selectMenu',
			'click button[data-button ^= "data-button"]': 'buttonClick',
			'blur input': 'inputDone',
			'keydown input': 'onEnter'
		},


		initialize: function(o, moleculeList) {
			this.moleculeList = moleculeList;
			this.model && this.listenTo(this.model, 'change', this.render);
			this.listenTo(this.moleculeList, 'all', this.showMoleculeList);
		},


		render: function() {
			this.$el.html(this.template(this.model.toJSON()));
			this.showMoleculeList();
			return this;
		},

		showMoleculeList: function() {
			var molecules = this.moleculeList.models;
			var ul = this.$('ul[target $= "Ref"]').html('');
			for (i = 0; i < molecules.length; i++){
				ul.each(function(){
					$(this).attr('type') === molecules[i].get('type') &&
					$(this).append('<li><a href = "#" val = "' + molecules[i].get('id') + '">' + molecules[i].get('id') + '</a></li>'); 
				})
			}
		},

		setModel: function (model) {
			this.stopListening();
			this.model = model;
			this.listenTo(this.model, 'change', this.render);
			this.listenTo(this.moleculeList, 'all', this.showMoleculeList);
			this.render();
		}
	});


	/*******************
	 * ControlView.js *
	 *******************/
	app.ControlView = app.BaseView.extend({

		// Cache the template function for a single item.
		template: _.template($('#control-property-template').html()),

		// The DOM events specific to an item.
		events: {
			'click button:nth-child(2)': 'chooseControl',
			'click ul[role ^= "menu"] a': 'selectMenu',
			'blur input': 'inputDone',
			'keydown input': 'onEnter'
		},

		initialize: function() {
			this.listenTo(this.model, 'change', this.render);
			this.listenTo(this.model, 'destroy', this.remove);
		},

		// Re-render the titles of the todo item.
		render: function() {
			this.$el.html(this.template(this.model.toJSON()));
			return this;
		},

		chooseControl: function(e) {
			e.preventDefault();
			var el = $(e.currentTarget);
			var attr = el.text().trim();
			var options = this.model.get('options');
			var newOptions = {};
			for (i in options)
				newOptions[i] = options[i];

			newOptions[attr] = !newOptions[attr];
			this.model.set('options', newOptions);
		},

		show: function() {
			this.$el.show();
		}
	});



	/***************
	 *  PTView.js  *
	 ***************/

	app.PTView = app.BaseView.extend({
		tagName: 'div',
		// Cache the template function for a single item.
		template: _.template($('#pt-template').html()),

		// The DOM events specific to an item.
		events: {
			'blur input': 'inputDone',
			'keydown input': 'onEnter',
			'click button[action ^= "delete"]': 'del',
            'click ul[role ^= "menu"] a': 'selectMenu'
		},

		initialize: function() {
			this.listenTo(this.model, 'change', this.render);
			this.listenTo(this.model, 'destroy', this.remove);
			this.render();
		},
		
		ptInputDone: function(e) {
			if (!this.inputDone.apply(this, [e])){
				this.model.destroy();
			}
			
		},
		
		render: function() {
			this.$el.html(this.template(this.model.toJSON()));
			return this;
		},

		show: function() {
			this.$el.show();
		},

		del: function() {
			this.model.destroy();
		}
	});


	/***************
	 *  PTView.js  *
	 ***************/
	app.ConditionsView = app.BaseView.extend({

		// Cache the template function for a single item.
		template: _.template($('#conditions-property-template').html()),

		// The DOM events specific to an item.
		events: {
			'click .group-head button': 'newPT',
			'click ul[role ^= "menu"] a': 'selectMenu',
            'blur input': 'inputDone',
			'keydown input': 'onEnter'
		},

		initialize: function(o, moleculeList) {
			this.moleculeList = moleculeList;
			this.pts = this.model.get('PTs');
			this.listenTo(this.model, 'change', this.render);
			this.listenTo(this.model, 'destroy', this.remove);
			this.listenTo(this.pts, 'add', this.addPT);
            this.listenTo(this.pts, 'reset', this.reset);
			this.listenTo(this.moleculeList, 'all', this.showMoleculeList);
		},

		render: function() {
			this.$el.html(this.template(this.model.toJSON()));
			this.showMoleculeList();
            var $el = this.$('.col-xs-12'),
                models = this.pts.models;
            var view, el, j = 0;
            for (var i = 0; i < models.length; i++){
	            if (models[i].get('P') || models[i].get('T')){
		            view = new app.PTView({model: models[i]});
                	$el.append(view.render().el);
	            }
            }
			return this;
		},

        reset: function(){
            this.render();
        },

		newPT: function(e) {
			var ptList = this.model.get('PTs');
			ptList.add(new app.PT());
		},

		addPT: function(model) {
//            var $left = this.$('.col-xs-6:first'),
//                $right = this.$('.col-xs-6:nth-child(2)'),
//                left_len = $left.children().length,
//                right_len = $right.children().length,
//                $el = left_len <= right_len ? $left : $right,
            var $el = this.$('.col-xs-12');
			var view = new app.PTView({model: model});
			$el.append(view.render().el);
			view.$el.find('input:first-child').focus();
		},

		showMoleculeList: function() {
			var molecules = this.moleculeList.models,
                $ul = this.$('ul[target $= "bathGas"]').html('');
			for (i = 0; i < molecules.length; i++){
				$ul.each(function(){
					$(this).attr('type') === molecules[i].get('type') &&
					$(this).append('<li><a href = "#" val = "' + molecules[i].get('id') + '">' + molecules[i].get('id') + '</a></li>');
				})
			}

            $ul = this.$('ul[target $= "initialPopulationRef"]').html('');
            for (i = 0; i < molecules.length; i++){
                if (molecules[i].get('type') != 'transitionState' && molecules[i].get('type') != 'bathGas'){
                    $ul.append('<li><a href = "#" val = "' + molecules[i].get('id') + '">' + molecules[i].get('id') + '</a></li>');
                }
			}
		},

		show: function() {
			this.$el.show();
		}
	});



	/***************
	 *  PTView.js  *
	 ***************/

	app.ModelParametersView = app.BaseView.extend({
		tagName: 'div',
		// Cache the template function for a single item.
		template: _.template($('#modelParameters-property-template').html()),

		// The DOM events specific to an item.
		events: {
			'blur input': 'inputDone',
			'keydown input': 'onEnter'
		},

		initialize: function() {
			this.listenTo(this.model, 'change', this.render);
			this.listenTo(this.model, 'destroy', this.remove);
		},

		render: function() {
			this.$el.html(this.template(this.model.toJSON()));
			return this;
		},

		show: function() {
			this.$el.show();
		}
	});

    $('#save').click(function(e){
        var postData = JSON.stringify(app.curData.toJSON());
        $.ajax({
            url: '/toxml/',
            async: true,
            type: 'POST',
            data: postData,
            success: function(data, type){
                var href = 'download/' + data;
                var iframe = document.getElementById('invisible');
                iframe.src = href;
           }
       })
    });

    $("#calculate").click(function(e){
        var postData = JSON.stringify(app.curData.toJSON());
        $.ajax({
            url: '/calculate/',
            async: true,
            type: 'POST',
            data: postData,
            success: function(data, type){
                var r = JSON.parse(data);
                if (r.error){
                    alert(r.error)
                } else {
                    var href = 'download/' + r.data;
                    var iframe = document.getElementById('invisible');
                    iframe.src = href;
                }

           }
       });
    });

    $('#openfile').on('change', onchange);

    function onchange(){
        var $this = $(this);
        if ($this.val()){
            $.ajaxFileUpload({
                url: '/openxml/',
                secureuri: false,
                fileElementId: 'openfile',
                dataType: 'json',
                success: function (data, status){
                    data = JSON.parse(data.replace('&gt;', '>'))
                    if(typeof(data.error) != 'undefined'){
                        alert(data.error);
                    } else {
                        app.curData.fromJSON(data);
                    }
                },
                error: function (data, status, e){
                    alert(e);
                }
            })
        }
        $this.replaceWith('<input class="btn btn-default" id="openfile" name="openfile" type="file" />');
        $('#openfile').on('change', onchange);
    }

    function downloadFile(fileName, content){
        var aLink = document.createElement('a');
        var blob = new Blob([content]);
        var evt = document.createEvent("HTMLEvents");
        evt.initEvent("click", false, true);
        aLink.download = fileName;
        aLink.href = URL.createObjectURL(blob);
        aLink.dispatchEvent(evt);
    }

    $('#renameModal button[action ^= "ok"]').click(function (e) {
        var $filename = $('#renameModal input[target ^= "filename"]');
        var $title = $('#renameModal input[target ^= "title"]');
        app.curData.set('filename', $filename.val());
        app.curData.set('title', $title.val());
        $('#renameModal').modal('hide');
    });

    $('#renameModal').on('show.bs.modal', function (e) {
        $('#renameModal input[target ^= "filename"]').val(app.curData.get('filename'));
        $('#renameModal input[target ^= "title"]').val(app.curData.get('title'));
    });

    $('#new').attr('href', location.href);
})(jQuery);

