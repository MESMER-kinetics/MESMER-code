var app = app || {};

app.BaseView = Backbone.View.extend({
	inputDone: function(e) {
		var $input = $(e.currentTarget);
		var value = $input.val().trim();
		var target = $input.attr('target');
		$input.val(value);
		if (value) {
			this.model.set(target, value);
			return value;
		}
		return false;
	},

	onEnter: function(e) {
		if (e.which === 13) {
			return this.inputDone(e);
		}
	},

	showMessage: function(message, type, el){
		el || (el = this.$el);
	},

	buttonClick: function(e) {
		var $el = $(e.currentTarget);
		var value = $el.attr('val');
		var target = $el.attr('target');
		this.model.set(target, value);
	},

	selectMenu: function(e) {
		e.preventDefault();
		var el = $(e.currentTarget);
		var key = el.parent().parent().attr('target');
		var value = el.attr('val');
		el.parent().parent().prev().dropdown('toggle');
		this.model.set(key, value);
		e.stopPropagation();
	}
})
